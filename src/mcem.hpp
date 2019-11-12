/** mccbn: large-scale inference on conjunctive Bayesian networks
 *
 *  This file is part of the mccbn package
 *
 * @author Susana Posada Céspedes
 * @email susana.posada@bsse.ethz.ch
 */

#ifndef MCEM_HPP
#define MCEM_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <random>
#include <vector>
#include <memory>
#include "rng_utils.hpp"

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
typedef Map<MatrixXd> MapMatd;
typedef Map<RowVectorXd> MapRowVecd;

struct Event {
  unsigned int event_id;
};

typedef boost::adjacency_list<
  boost::hash_setS, boost::vecS, boost::bidirectionalS, Event> Poset;
typedef boost::graph_traits<Poset>::vertices_size_type vertices_size_type;
typedef boost::graph_traits<Poset>::vertex_descriptor Node;
typedef std::pair<vertices_size_type, vertices_size_type> Edge;
typedef std::vector<Node> node_container;
typedef std::vector<Edge> edge_container;


class Context {
public:
  typedef ::rng_type rng_type;
  typedef boost::ptr_vector<rng_type> rng_vector_type;

  rng_type rng;

  Context(int seed, bool verbose=false): rng(seed), _verbose(verbose) {}

  inline bool get_verbose() const {
    return _verbose;
  }

  std::unique_ptr<rng_vector_type> get_auxiliary_rngs(int num_rngs) {
    std::unique_ptr<rng_vector_type> rngs(new rng_vector_type());
    rngs->reserve(num_rngs);
    for (int t = 0; t < num_rngs; ++t)
      rngs->push_back(new rng_type(rng()));
    return std::move(rngs);
  }

protected:
  bool _verbose;
};

class Model {
public:
  Poset poset;                // Adjacency list encoding cover releations
  node_container topo_path;   // A topological ordering of the nodes
  bool cycle;                 // Boolean variable indicating whether the poset contains cycles
  bool reduction_flag;        // Boolean variable indicating whether the poset is a transitively reduced graph

  /* Default constructor */
  Model(float lambda_s=1.0, bool cycle=false, bool reduction=false,
        bool update_children=true, bool update_node_idx=true) : cycle(cycle),
        reduction_flag(reduction), _lambda_s(lambda_s),
        _update_children(update_children), _update_node_idx(update_node_idx),
        _size(0) {}

  /* Parametrized constructor */
  Model(unsigned int p, float lambda_s=1.0, bool cycle=false,
        bool reduction=false, bool update_children=true,
        bool update_node_idx=true) : poset(p), cycle(cycle),
        reduction_flag(reduction), _lambda(p), _lambda_s(lambda_s),
        _children(p), _update_children(update_children), _node_idx(p),
        _update_node_idx(update_node_idx), _size(p) {
    topo_path.reserve(p);
  }

  /* Constructor using the edge iterator constructor */
  Model(const edge_container& edge_list, unsigned int p, float lambda_s=1.0,
         bool cycle=false, bool reduction=false, bool update_children=true) :
    poset(edge_list.begin(), edge_list.end(), p), cycle(cycle),
    reduction_flag(reduction), _lambda(p), _lambda_s(lambda_s),
    _children(p), _update_children(update_children), _node_idx(p), _size(p) {
    topo_path.reserve(p);
    for (unsigned int i = 0; i < p; ++i) {
      poset[i].event_id = i;
      _node_idx[i] = i;
      _update_node_idx = false;
    }
  }

  /* Copy constructor */
  Model(const Model& m) : poset(m.poset), topo_path(m.topo_path),
  cycle(m.cycle), reduction_flag(m.reduction_flag), _lambda(m.get_lambda()),
  _lambda_s(m.get_lambda_s()), _epsilon(m.get_epsilon()),
  _llhood(m.get_llhood()), _children(m.get_children()),
  _update_children(m.get_update_children()), _node_idx(m.get_node_idx()),
  _update_node_idx(m.get_update_node_idx()), _size(m.size()) {}

  inline vertices_size_type size() const;

  void set_lambda(const Eigen::Ref<const VectorXd>& lambda);

  void set_lambda(const Eigen::Ref<const VectorXd>& lambda,
                  const float max_lambda);

  void set_epsilon(const double eps);

  void set_llhood(const double llhood);

  void set_children();

  inline VectorXd get_lambda() const;

  inline double get_lambda(const unsigned int idx) const;

  inline float get_lambda_s() const;

  inline double get_epsilon() const;

  inline double get_llhood() const;

  std::vector<node_container> get_direct_predecessors() const;

  inline const std::vector< std::unordered_set<Node> >& get_children() const;

  inline const std::vector<int>& get_node_idx() const;

  inline const int get_node_idx(const unsigned int i) const;

  inline bool get_update_children() const;

  inline bool get_update_node_idx() const;

  void has_cycles();

  void topological_sort();

  void update_children(const Node& u);

  std::vector<node_container> get_direct_successors(node_container& topo_order);

  void transitive_reduction_dag();

  bool add_edge(const Node& u, const Node& v);

  void remove_edge(const Node& u, const Node& v);

  void remove_redundant_edges(const Node& u, const Node& v);

  bool add_relation(const Node& u, const Node& v);

  void remove_relation(const Node& u, const Node& v);

  void swap_node(const Node& u, const Node& v);

  void update_node_idx();

  template <typename PropertyMap>
  void print_cover_relations(PropertyMap name);

  void print_cover_relations();

  void clear();

protected:
  VectorXd _lambda;
  float _lambda_s;
  double _epsilon;
  double _llhood;
  std::vector< std::unordered_set<Node> > _children;
  bool _update_children;
  std::vector<int> _node_idx;
  bool _update_node_idx;
  vertices_size_type _size;
};

class DataImportanceSampling {
public:
  VectorXd w;
  VectorXi dist;
  MatrixXd Tdiff;

  /*Parametrized Constructor */
  DataImportanceSampling(unsigned int L, unsigned int p) : w(L), dist(L),
    Tdiff(L, p) {}

};

/* Class containing customisable options for the EM algorithm */
class ControlEM {
public:
  unsigned int max_iter;
  unsigned int update_step_size; // increase L every 'update_step_size' to reach a desirable 'tol'
  double tol;                    // convergence tolerance
  float max_lambda;
  unsigned int neighborhood_dist;

  ControlEM(unsigned int max_iter=100, unsigned int update_step_size=20,
            double tol=0.001, float max_lambda=1e6,
            unsigned int neighborhood_dist=1) :
    max_iter(max_iter), update_step_size(update_step_size), tol(tol),
    max_lambda(max_lambda), neighborhood_dist(neighborhood_dist) {}
};

vertices_size_type Model::size() const {
  return _size;
}

VectorXd Model::get_lambda() const {
  return _lambda;
}

double Model::get_lambda(const unsigned int idx) const {
  return _lambda[idx];
}

float Model::get_lambda_s() const {
  return _lambda_s;
}

double Model::get_epsilon() const {
  return _epsilon;
}

double Model::get_llhood() const {
  return _llhood;
}

const std::vector< std::unordered_set<Node> >& Model::get_children() const {
  return _children;
}

const std::vector<int>& Model::get_node_idx() const {
  return _node_idx;
}

const int Model::get_node_idx(const unsigned int i) const {
  return _node_idx[i];
}


bool Model::get_update_children() const {
  return _update_children;
}

bool Model::get_update_node_idx() const {
  return _update_node_idx;
}

DataImportanceSampling importance_weight(
    const RowVectorXb& genotype, unsigned int L, const Model& model,
    const double time, const std::string& sampling,
    const VectorXd& scale_cumulative, const VectorXi& dist_pool,
    const MatrixXd& Tdiff_pool, const unsigned int neighborhood_dist,
    Context::rng_type& rng, const bool sampling_times_available=false);

VectorXi hamming_dist_mat(const MatrixXb &x, const RowVectorXb &y);

double complete_log_likelihood(
    const VectorXd &lambda, const double eps, const MatrixXd &Tdiff,
    const VectorXd &dist, const float W, const bool internal=true);

double MCEM_hcbn(
    Model& model, const MatrixXb& obs, const VectorXd& times,
    const RowVectorXd& weights, const unsigned int L,
    const std::string& sampling, const ControlEM& control_EM,
    const bool sampling_times_available, const unsigned int thrds,
    Context& ctx);

unsigned int n_choose_k(unsigned int n, unsigned int k);

void neighbors(const unsigned int p, unsigned int k,
               Eigen::Ref<MatrixXb> output_mat);

void coin_tossing(MatrixXb& output_mat, const double eps,
                  Context::rng_type& rng);

edge_container adjacency_mat2list(const MatrixXi& poset);

MatrixXi adjacency_list2mat(const Model& model);

bool is_compatible(const RowVectorXb& genotype, const Model& model);

int num_compatible_observations(const MatrixXb& obs, const Model& poset,
                                const unsigned int N);

int num_incompatible_events(const MatrixXb& genotype, const Model& poset);

std::vector<int> rdiscrete_std(const unsigned int N, const VectorXd& weights,
                               Context::rng_type& rng);

void handle_exceptions();

Rcpp::NumericMatrix generateMutationTimes(
    const Rcpp::NumericVector& lambda, const Rcpp::IntegerMatrix& poset,
    Rcpp::IntegerMatrix& O, const std::vector<double>& times,
    const Rcpp::IntegerVector& topo_path, int nrOfSamples, float lambda_s,
    bool sampling_time_available);

#endif
