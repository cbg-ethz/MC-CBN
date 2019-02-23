#ifndef MCEM_H
#define MCEM_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <vector>

using Eigen::Map;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::RowVectorXd;
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;
typedef Eigen::Matrix<bool, 1, Eigen::Dynamic> RowVectorXb;
typedef Map<VectorXd> MapVecd;
typedef Map<VectorXi> MapVeci;
typedef Map<MatrixXi> MapMati;
typedef Map<RowVectorXd> MapRowVecd;

class Context {
public:
  typedef std::mt19937 rng_type;
  typedef std::vector<rng_type> rng_vector_type;

  rng_type rng;

  Context(int seed, bool verbose=false): _verbose(verbose) {
    _set_seed(rng, seed);
  }

  bool get_verbose() const {
    return _verbose;
  }

  rng_type get_auxiliary_rng() {
    rng_type aux_rng;
    _set_seed(aux_rng, rng());
    return aux_rng;
  }

  rng_vector_type get_auxiliary_rngs(int num_rngs) {
    std::vector<rng_type> rngs;
    rngs.reserve(num_rngs);
    for (int t = 0; t < num_rngs; ++t)
      rngs.push_back(get_auxiliary_rng());
    return rngs;
  }

protected:
  bool _verbose;

  static void _set_seed(rng_type& rng, int seed) {
    std::ranlux24_base seeder_rng = std::ranlux24_base(seed);

    // adapted from https://stackoverflow.com/a/15509942
    std::array<int, std::mt19937::state_size> seed_data;
    std::generate_n(seed_data.data(), seed_data.size(), std::ref(seeder_rng));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));

    rng.seed(seq);
  }
};

class Model {
public:
  MatrixXi poset;                      // Matrix encoding cover relations
  VectorXi N_pa;                       // Number of parents for each event/mutation
  std::vector< std::vector<int> > pa;  // List of parents per event
  VectorXi N_pa_closure;               // Number of parents for a transitively closed graph
  std::vector< std::vector<int> > pa_closure; // List of parents, grandparents, great-grandparents... per event
  std::vector< std::vector<int> > ch_closure; // List of all the childreen per event

  // Parametrized Constructors
  Model(unsigned int p, float lambda_s=1.0) : poset(p, p), N_pa(p), pa(p),
    _lambda(p), _lambda_s(lambda_s) {}

  Model(unsigned int p, bool asa, float lambda_s=1.0) :
    poset(p, p), N_pa(p), pa(p), N_pa_closure(p), pa_closure(p), ch_closure(p),
    _lambda(p), _lambda_s(lambda_s) {}

  unsigned int size() const;

  void fill_poset(const Eigen::Ref<const MatrixXi>& P);

  void set_lambda(const Eigen::Ref<const VectorXd>& lambda);

  void set_lambda(const Eigen::Ref<const VectorXd>& lambda,
                  const float max_lambda);

  void set_epsilon(const double eps);

  void set_llhood(const double llhood);

  VectorXd get_lambda() const;

  double get_lambda(const unsigned int idx) const;

  float get_lambda_s() const;

  double get_epsilon() const;

  double get_llhood() const;

  void grandparents(unsigned int current, unsigned int parent,
                    std::vector<bool>& visit);

  void grandchildren(unsigned int current, unsigned int child,
                     std::vector<bool>& visit);

  void parents();

  void parents(const Eigen::Ref<const VectorXi>& topological_path);

  void children(const Eigen::Ref<const VectorXi>& topological_path);

protected:
  VectorXd _lambda;
  float _lambda_s;
  double _epsilon;
  double _llhood;
};

class DataImportanceSampling {
public:
  VectorXd w;
  VectorXi dist;
  MatrixXd Tdiff;

  //Parametrized Constructor
  DataImportanceSampling(unsigned int L, unsigned int p) : w(L), dist(L),
    Tdiff(L, p) {}

};

DataImportanceSampling importance_weight(
    const RowVectorXb &genotype, const unsigned int L, const Model& model,
    const VectorXi &topo_path, const double time, const std::string& sampling,
    const unsigned int version, const float perturb_prob,
    const VectorXi &dist_pool, const MatrixXd &Tdiff_pool,
    Context::rng_type& rng, const bool sampling_times_available);

VectorXi hamming_dist_mat(const MatrixXb &x, const RowVectorXb &y);

double complete_log_likelihood(
    const VectorXd &lambda, const double eps, const MatrixXd &Tdiff,
    const VectorXd &dist, const float W);

double MCEM_hcbn(
    Model& model, const MatrixXb &obs, const VectorXd& times,
    const VectorXi& topo_path, const RowVectorXd& weights,
    const unsigned int L, const std::string& sampling,
    const unsigned int version, const float perturb_prob,
    const unsigned int max_iter, const float burn_in, const float max_lambda,
    const bool sampling_times_available,
    const unsigned int thrds, Context& ctx);

Rcpp::NumericMatrix generateMutationTimes(
    const Rcpp::NumericVector &lambda, const Rcpp::IntegerMatrix &poset,
    Rcpp::IntegerMatrix &O, const std::vector<double> &times,
    const Rcpp::IntegerVector &topo_path, int nrOfSamples, float lambda_s,
    bool sampling_time_available);

#endif
