/** mccbn: large-scale inference on conjunctive Bayesian networks
 *  MCEM for the hidden conjuctive Bayesian network model
 *
 * @author Susana Posada Céspedes
 * @email susana.posada@bsse.ethz.ch
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "mcem.hpp"
#include "add_remove.hpp"
#include "not_acyclic_exception.hpp"
#include <boost/graph/graph_traits.hpp>
#include <random>
#include <vector>

// #include "debugging_helper.hpp"

#ifdef _OPENMP
  #include <omp.h>
#endif

class Initializer {
public:
  Initializer() {
    Eigen::initParallel();
  }
} initializer;

void handle_exceptions() {
  try {
    throw;
  } catch (const std::exception& ex) {
    // NOTE: reference to 'exception' is ambiguous, 'std::' required
    std::string msg = std::string("c++ exception: ") + ex.what();
    ::Rf_error(msg.c_str());
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}

//' @noRd
//' @param N number of samples to be drawn
std::vector<int> rdiscrete_std(const unsigned int N, const VectorXd& weights,
                               Context::rng_type& rng) {
  std::discrete_distribution<int> distribution(weights.data(),
                                               weights.data() + weights.size());
  std::vector<int> ret(N);

  for (unsigned int i = 0; i < N; ++i)
    ret[i] = distribution(rng);

  return ret;
}

VectorXd log_bernoulli_process(const VectorXd& dist, const double eps,
                               const unsigned int p) {

  const unsigned int L = dist.size();
  VectorXd log_prob = VectorXd::Zero(L);

  if (eps == 0) {
    /* NOTE: If all observations are compatible with poset (which can happen
     * because noisy observations can be compatible), then eps can be 0, as
     * well as dist.
     */
    for (unsigned int i = 0; i < L; ++i) {
      if (dist[i] != 0) {
        double log_eps = std::log(eps + DBL_EPSILON);
        double log1m_eps = std::log1p(- eps - DBL_EPSILON);
        log_prob[i] = log_eps * dist[i] + log1m_eps * (p - dist[i]);
      } else {
        log_prob[i] = 0;
      }
    }
  } else {
    double log_eps = std::log(eps);
    double log1m_eps = std::log1p(-eps);
    log_prob = (log_eps * dist).array() + log1m_eps * (p - dist.array());
  }
  return log_prob;
}

double log_bernoulli_process(const unsigned int dist, const double eps,
                             const unsigned int p) {
  double log_prob = 0.0;

  if (eps == 0) {
    /* NOTE: If all observations are compatible with poset (which can happen
    * because noisy observations can be compatible), then eps can be 0, as
    * well as dist.
    */
    if (dist != 0)
      log_prob = std::log(eps + DBL_EPSILON) * dist +
        std::log1p(- eps - DBL_EPSILON) * (p - dist);
  } else {
    log_prob = std::log(eps) * dist + std::log1p(-eps) * (p - dist);
  }
  return log_prob;
}

//' Compute complete-data log-likelihood or (equivalently) hidden log-likelihood
double complete_log_likelihood(const VectorXd& lambda, const double eps,
                               const MatrixXd& Tdiff, const VectorXd& dist,
                               const float W, const bool internal) {

  const unsigned int p = lambda.size();
  const unsigned int N = dist.size();
  double llhood;

  switch(internal) {
  case(0):
    llhood = N * lambda.array().log().sum() - (Tdiff * lambda).sum();
    break;
  case(1):
    llhood = N * (lambda.array().log().sum() - p);
  }

  if (eps == 0) {
    for (unsigned int i = 0; i < N; ++i) {
      if (dist(i) != 0)
        llhood += std::log(eps + DBL_EPSILON) * dist(i) +
          std::log(1 - eps - DBL_EPSILON) * (p - dist(i));
    }
  } else {
    llhood += std::log(eps) * dist.sum() + std::log1p(- eps) * (p - dist.array()).sum();
  }
  return llhood;
}

//' Generate observations from a given poset and given rates
//'
//' @noRd
//' @param N number of samples
//' @return returns matrix containing observations
MatrixXb sample_genotypes(
    const unsigned int N, const Model& model, MatrixXd& T_events,
    VectorXd& T_sampling, Context::rng_type& rng,
    const bool sampling_times_available=false) {

  /* Initialization and instantiation of variables */
  const vertices_size_type p = model.size();  // Number of mutations / events
  MatrixXd T_events_sum;
  MatrixXb obs;
  T_events_sum.setZero(N, p);
  obs.setConstant(N, p, false);
  auto id = boost::get(&Event::event_id, model.poset);

  /* Generate occurence times T_events_{j} ~ Exp(lambda_{j}) */
  for (unsigned int j = 0; j < p; ++j)
    T_events.col(j) = rexp(N, model.get_lambda(j), rng);


  /* Use sampling times when available */
  if (!sampling_times_available)
    T_sampling = rexp(N, model.get_lambda_s(), rng);

  VectorXd T_max = VectorXd::Zero(N);
  /* Loop through nodes in topological order */
  for (node_container::const_reverse_iterator v = model.topo_path.rbegin();
       v != model.topo_path.rend(); ++v) {
    T_max = VectorXd::Zero(N);
    /* Loop through parents of node v */
    boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
    for (boost::tie(in_begin, in_end) = boost::in_edges(*v, model.poset);
         in_begin != in_end; ++in_begin) {
      Node u = boost::get(id, source(*in_begin, model.poset));
      T_max = (T_events_sum.col(u).array() > T_max.array()).select(T_events_sum.col(u), T_max);
    }
    T_events_sum.col(model.poset[*v].event_id) = T_events.col(model.poset[*v].event_id) + T_max;
    obs.col(model.poset[*v].event_id) =
      (T_events_sum.col(model.poset[*v].event_id).array() <= T_sampling.array()).select(true, obs.col(model.poset[*v].event_id));
  }
  return obs;
}

//' Generate observation times from a given poset and given rates
//'
//' @noRd
//' @param N number of samples
//' @return returns matrix containing occurrence times
MatrixXd sample_times(
    const unsigned int N, const Model& model, MatrixXd& T_events,
    Context::rng_type& rng) {

  /* Initialization and instantiation of variables */
  const vertices_size_type p = model.size();  // Number of mutations / events
  MatrixXd T_events_sum;
  T_events_sum.setZero(N, p);
  auto id = boost::get(&Event::event_id, model.poset);

  /* Generate occurence times T_events_{j} ~ Exp(lambda_{j}) */
  for (unsigned int j = 0; j < p; ++j)
    T_events.col(j) = rexp(N, model.get_lambda(j), rng);

  VectorXd T_max = VectorXd::Zero(N);
  /* Loop through nodes in topological order */
  for (node_container::const_reverse_iterator v = model.topo_path.rbegin();
       v != model.topo_path.rend(); ++v) {
    T_max = VectorXd::Zero(N);
    /* Loop through parents of node v */
    boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
    for (boost::tie(in_begin, in_end) = boost::in_edges(*v, model.poset);
         in_begin != in_end; ++in_begin) {
      Node u = boost::get(id, source(*in_begin, model.poset));
      T_max = (T_events_sum.col(u).array() > T_max.array()).select(T_events_sum.col(u), T_max);
    }
    T_events_sum.col(model.poset[*v].event_id) = T_events.col(model.poset[*v].event_id) + T_max;
  }
  return T_events_sum;
}

MatrixXb generate_genotypes(
    const MatrixXd& T_events_sum, const Model& model,
    VectorXd& T_sampling, Context::rng_type& rng,
    const bool sampling_times_available=false) {

  /* Initialization and instantiation of variables */
  const unsigned int N = T_events_sum.rows(); // Number of genotypes to be drawn
  const vertices_size_type p = model.size();  // Number of mutations / events
  MatrixXb obs;
  obs.setConstant(N, p, false);

  if (!sampling_times_available)
    T_sampling = rexp(N, model.get_lambda_s(), rng);

  /* Loop through nodes in topological order */
  for (node_container::const_reverse_iterator v = model.topo_path.rbegin();
       v != model.topo_path.rend(); ++v)
    obs.col(model.poset[*v].event_id) =
      (T_events_sum.col(model.poset[*v].event_id).array() <= T_sampling.array()).select(true, obs.col(model.poset[*v].event_id));

  return obs;
}

//' Compute observed log-likelihood
double obs_log_likelihood(
    const MatrixXb& obs, const MatrixXi& poset, const VectorXd& lambda,
    const double eps, const RowVectorXd& weights, const VectorXd& times,
    const unsigned int L, const std::string& sampling,
    const unsigned int neighborhood_dist, Context& ctx, const float lambda_s=1.0,
    const bool sampling_times_available=false, const unsigned int thrds=1) {

  const auto p = poset.rows(); // Number of mutations / events
  const auto N = obs.rows();   // Number of observations / genotypes
  double llhood = 0.0;

  edge_container edge_list = adjacency_mat2list(poset);
  Model model(edge_list, p, lambda_s);
  model.set_lambda(lambda);
  model.set_epsilon(eps);
  model.has_cycles();
  if (model.cycle) {
    throw not_acyclic_exception();
  } else {
    model.topological_sort();

    unsigned int K = 0;
    VectorXd scale_cumulative;
    MatrixXd Tdiff_pool;
    MatrixXd T_pool;
    if (sampling == "add-remove") {
      scale_cumulative.resize(p);
      scale_cumulative = scale_path_to_mutation(model);
      if (model.get_update_node_idx())
        model.update_node_idx();
    } else if (sampling == "pool") {
      K = p * L;
      Tdiff_pool.resize(K, p);
      T_pool.resize(K, p);
      T_pool = sample_times(K, model, Tdiff_pool, ctx.rng);
    }

    #ifdef _OPENMP
    omp_set_num_threads(thrds);
    #endif
    auto rngs = ctx.get_auxiliary_rngs(thrds);

    #pragma omp parallel for reduction(+:llhood) schedule(static)
    for (unsigned int i = 0; i < N; ++i) {
      VectorXi d_pool;
      if (sampling == "pool") {
        VectorXd T_sampling(K);
        if (sampling_times_available)
          T_sampling.setConstant(times[i]);

        MatrixXb genotype_pool = 
          generate_genotypes(T_pool, model, T_sampling, (*rngs)[omp_get_thread_num()],
                             sampling_times_available);
        d_pool = hamming_dist_mat(genotype_pool, obs.row(i));
      }
      DataImportanceSampling importance_sampling = importance_weight(
        obs.row(i), L, model, times[i], sampling, scale_cumulative, d_pool,
        Tdiff_pool, neighborhood_dist, (*rngs)[omp_get_thread_num()],
        sampling_times_available);

      double aux = weights(i) * importance_sampling.w.sum();
      if (aux > 0) {
        int L_eff = L;
        if (sampling == "backward")
          L_eff = (importance_sampling.w.array() > 0).count();
        llhood += std::log(aux / L_eff);
      }
    }
  }
  return llhood;
}

//' Compute Hamming distance between two vectors
//'
//' @noRd
//' @return returns Hamming distance
int hamming_dist(const VectorXi& x, const VectorXi& y) {
  return (x - y).array().abs().sum();
}

//' Compute Hamming distance between a matrix and a vector row-wise
//'
//' @noRd
//' @return returns a vector containing the Hamming distance
VectorXi hamming_dist_mat(const MatrixXb& x, const RowVectorXb& y) {
  const int N = x.rows();
  return (x.array() != y.replicate(N, 1).array()).rowwise().count().cast<int>();
  // return (x.rowwise() - y).array().abs().rowwise().sum();
}

//' Compute importance weights and (expected) sufficient statistics by
//' importance sampling
//'
//' @noRd
//' @param genotype a p-dimesional vector corresponding to an observed
//' genotype - e.g., mutated (1) and non-mutated (0) genes
//' @param L number of samples
//' @param model
//' @param time sampling time
//' @param sampling variable indicating which proposal to use
//' @return returns importance weights and (expected) sufficient statistics
DataImportanceSampling importance_weight(
    const RowVectorXb& genotype, unsigned int L, const Model& model,
    const double time, const std::string& sampling,
    const VectorXd& scale_cumulative, const VectorXi& dist_pool,
    const MatrixXd& Tdiff_pool, const unsigned int neighborhood_dist,
    Context::rng_type& rng, const bool sampling_times_available) {

  /* Initialization and instantiation of variables */
  const vertices_size_type p = model.size(); // Number of mutations / events
  unsigned int reps = 0;
  unsigned int L_aux = 1;
  std::vector<int> nrows(neighborhood_dist);
  unsigned int nrows_sum = 0;
  if (sampling == "backward") {
    /* Using the leading and the first k order terms in X */
    for (unsigned int i = 0; i <= neighborhood_dist; ++i) {
      nrows[i] = n_choose_k(p, i);
      nrows_sum += nrows[i];
    }
    reps = L / nrows_sum;
    L = reps * nrows_sum;
  } else if (sampling == "bernoulli") {
    reps = 1;
    L_aux = std::max(L / reps, L_aux);
    L = reps * L_aux;
  }
  DataImportanceSampling importance_sampling(L, p);

  if (sampling == "forward") {
    /* Generate L samples from poset with parameters lambda and lambda_s.
     * In particular, epsilon is zero (default value) - because the idea is to
     * generate samples of X (true genotype)
     */
    MatrixXb samples(L, p);
    VectorXd T_sampling(L);
    if (sampling_times_available)
      T_sampling.setConstant(time);

    samples = sample_genotypes(L, model, importance_sampling.Tdiff, T_sampling,
                               rng, sampling_times_available);
    importance_sampling.dist = hamming_dist_mat(samples, genotype);
    // importance_sampling.w = pow(model.get_epsilon(), importance_sampling.dist.array()) *
    //   pow(1 - model.get_epsilon(), p - importance_sampling.dist.array());
    VectorXd d = importance_sampling.dist.cast<double>();
    importance_sampling.w = pow(model.get_epsilon(), d.array()) *
      pow(1 - model.get_epsilon(), p - d.array());
  } else if (sampling == "add-remove") {
    MatrixXb samples(L, p);
    VectorXd log_prob_Y_X(L);
    VectorXd log_prob_X(L);
    VectorXd log_proposal(L);

    /* Pick a move: add or remove with equal probability
     * Remove: choose event x with prob. ~ sum_{j \in max_path to x}(1/lambda_j)
     * (not applicable if genotype corresponds to the wild type)
     * Add: choose event x with an inverse probability (not applicable if the
     * genotype corresponds to the resistant genotype)
     * 'Stand-still': remain unperturbed (applicable if genotype is compatible
     * with the poset)
     */
    bool compatible = is_compatible(genotype, model);
    unsigned int mutations = genotype.count();
    bool wild_type = (mutations == 0);
    bool resistant_type = (mutations == p);
    std::vector<int> move(L);
    if (compatible) {
      if (wild_type) {
        /* possible moves: add (move 0) or stand-still (move 2) */
        Eigen::Vector3d weights(0.5, 0.0, 0.5);
        move = rdiscrete_std(L, weights, rng);
        log_proposal.setConstant(std::log(0.5));
      } else if (resistant_type) {
        /* possible moves: remove (move 1) or stand-still (move 2) */
        Eigen::Vector3d weights(0.0, 0.5, 0.5);
        move = rdiscrete_std(L, weights, rng);
        log_proposal.setConstant(std::log(0.5));
      } else {
        Eigen::Vector3d weights = Eigen::Vector3d::Ones();
        move = rdiscrete_std(L, weights, rng);
        log_proposal.setConstant(-std::log(3));
      }
    } else {
      /* possible moves: add (move 0) or remove (move 1) */
      Eigen::Vector2d weights = Eigen::Vector2d::Constant(0.5);
      move = rdiscrete_std(L, weights, rng);
      log_proposal.setConstant(std::log(0.5));
    }
    double q_choice;
    VectorXd remove_weight = genotype.select(scale_cumulative.transpose(), 0);
    VectorXd add_weight = scale_cumulative.array().inverse();
    add_weight = genotype.select(0, add_weight.transpose());
    std::vector<int> idx_remove = rdiscrete_std(L, remove_weight, rng);
    std::vector<int> idx_add = rdiscrete_std(L, add_weight, rng);

    for (unsigned int l = 0; l < L; ++l) {
      samples.row(l) = draw_sample(genotype, model, move[l], remove_weight,
                                   add_weight, q_choice, idx_remove[l],
                                   idx_add[l], compatible);
      log_proposal[l] += std::log(q_choice);
    }
    importance_sampling.dist = hamming_dist_mat(samples, genotype);

    VectorXd T_sampling(L);
    if (sampling_times_available)
      T_sampling.setConstant(time);

    /* Generate mutation times based on samples */
    importance_sampling.Tdiff =
      generate_mutation_times(samples, model, log_proposal, T_sampling, rng,
                              sampling_times_available);
    log_prob_Y_X =
      log_bernoulli_process(importance_sampling.dist.cast<double>(),
                            model.get_epsilon(), p);
    log_prob_X = cbn_density_log(importance_sampling.Tdiff, model.get_lambda());

    importance_sampling.w =
      (log_prob_Y_X + log_prob_X - log_proposal).array().exp();
  } else if (sampling == "pool") {
    // VectorXd q_prob = pow(model.get_epsilon(), dist_pool.array()) *
    //   pow(1 - model.get_epsilon(), p - dist_pool.array());
    VectorXd d_pool = dist_pool.cast<double>();
    VectorXd q_prob = pow(model.get_epsilon(), d_pool.array()) *
      pow(1 - model.get_epsilon(), p - d_pool.array());
    /* In the unlikely event that q_prob is 0 for all samples, default to
     * random sampling
     */
    bool random = false;
    if (q_prob.sum() == 0) {
      q_prob.setConstant(1);
      random = true;
    }
    double q_prob_sum = q_prob.sum();
    q_prob /= q_prob_sum;

    /* Draw L samples with replacement and with weights q_prob */
    std::vector<int> idxs_sample = rdiscrete_std(L, q_prob, rng);
    int idx;
    for (unsigned int l = 0; l < L; ++l) {
      idx = idxs_sample[l];
      importance_sampling.dist(l) = dist_pool(idx);
      importance_sampling.Tdiff.row(l) = Tdiff_pool.row(idx);
    }

    if (random)
      importance_sampling.w =
        log_bernoulli_process(importance_sampling.dist.cast<double>(),
                              model.get_epsilon(), p).array().exp();
    else
      importance_sampling.w.setConstant(q_prob_sum / dist_pool.size());
  } else if (sampling == "backward") {
    VectorXd log_prob_Y_X(L);
    VectorXd log_prob_X(L);
    VectorXd log_proposal = VectorXd::Zero(L);

    MatrixXb samples = genotype.replicate(nrows_sum, 1);
    VectorXi dist = VectorXi::Zero(nrows_sum);
    VectorXd log_prob_Y_X_aux = VectorXd::Zero(nrows_sum);
    log_prob_Y_X_aux[0] =
      log_bernoulli_process(0.0, model.get_epsilon(), p);
    nrows_sum = 1;
    for (unsigned int i = 1; i <= neighborhood_dist; ++i) {
      neighbors(p, i, samples.block(nrows_sum, 0, nrows[i], p));
      dist.segment(nrows_sum, nrows[i]).setConstant(i);
      log_prob_Y_X_aux.segment(nrows_sum, nrows[i]).setConstant(
          log_bernoulli_process(i, model.get_epsilon(), p));
      nrows_sum += nrows[i];
    }

    MatrixXb samples_rep = samples.replicate(reps, 1);
    importance_sampling.dist = dist.replicate(reps, 1);
    log_prob_Y_X = log_prob_Y_X_aux.replicate(reps, 1);

    VectorXd T_sampling(L);
    if (sampling_times_available)
      T_sampling.setConstant(time);

    /* Generate mutation times based on samples */
    importance_sampling.Tdiff =
      generate_mutation_times(samples_rep, model, log_proposal, T_sampling, rng,
                              sampling_times_available);
    /* Account for incompatible samples */
    Eigen::Matrix<bool, Eigen::Dynamic, 1> incompatible_samples = (importance_sampling.Tdiff.array() < 0.0).rowwise().any();
    log_proposal = log_proposal.array() - std::log(nrows_sum - incompatible_samples.count() / reps);

    log_prob_X = cbn_density_log(importance_sampling.Tdiff, model.get_lambda());

    importance_sampling.w =
      (log_prob_Y_X + log_prob_X - log_proposal).array().exp();
    /* Downweight samples that are not feasible / incompatible with current poset */
    importance_sampling.w = incompatible_samples.select(0, importance_sampling.w);
  } else if (sampling == "bernoulli") {
    VectorXd log_prob_Y_X(L);
    VectorXd log_prob_X(L);
    VectorXd log_proposal = VectorXd::Zero(L);

    MatrixXb samples = genotype.replicate(L_aux, 1);
    coin_tossing(samples, model.get_epsilon(), rng);

    VectorXi dist = VectorXi::Zero(L_aux);
    VectorXd log_prob_Y_X_aux = VectorXd::Zero(L_aux);
    dist = hamming_dist_mat(samples, genotype);

    MatrixXb samples_rep = samples.replicate(reps, 1);
    importance_sampling.dist = dist.replicate(reps, 1);

    VectorXd T_sampling(L);
    if (sampling_times_available)
      T_sampling.setConstant(time);

    /* Generate mutation times based on samples */
    importance_sampling.Tdiff =
      generate_mutation_times(samples_rep, model, log_proposal, T_sampling, rng,
                              sampling_times_available);

    log_prob_X = cbn_density_log(importance_sampling.Tdiff, model.get_lambda());

    importance_sampling.w = (log_prob_X - log_proposal).array().exp();
    /* Account for incompatible samples and downweight them */
    Eigen::Matrix<bool, Eigen::Dynamic, 1> incompatible_samples = (importance_sampling.Tdiff.array() < 0.0).rowwise().any();
    importance_sampling.w = incompatible_samples.select(0, importance_sampling.w);
  }

  return importance_sampling;
}

//' Compute importance weights and sufficient statistics by sampling
//'
//' @noRd
double MCEM_hcbn(
    Model& model, const MatrixXb& obs, const VectorXd& times,
    const RowVectorXd& weights, unsigned int L, const std::string& sampling,
    const ControlEM& control_EM, const bool sampling_times_available,
    const unsigned int thrds, Context& ctx) {

  // Initialization and instantiation of variables
  const vertices_size_type p = model.size(); // Number of mutations / events
  const unsigned int N = obs.rows();         // Number of observations / genotypes
  unsigned int N_eff = 0;
  unsigned int K = 0;
  unsigned int update_step_size = control_EM.update_step_size;
  VectorXd avg_lambda = VectorXd::Zero(p);
  VectorXd avg_lambda_current = VectorXd::Zero(p);
  double avg_eps = 0.0, avg_llhood = 0.0, obs_llhood = 0.0;
  double avg_eps_current = 0.0;
  bool tol_comparison = true;
  double expected_dist = 0.0;
  MatrixXd expected_Tdiff = MatrixXd::Zero(N, p);
  VectorXd Tdiff_colsum(p);
  MatrixXd Tdiff_pool;
  VectorXd scale_cumulative;

  if (sampling == "add-remove") {
    scale_cumulative.resize(p);
    if (model.get_update_node_idx())
      model.update_node_idx();
  } else if (sampling == "pool") {
    // if (p < 19) {
    //   K = std::max((unsigned int) pow(2, p + 1), 2 * L);
    // } else {
    //   // K = 100000000 / (sizeof(double) * p); //NOTE:empirically 100000 limit on runtime
    //   K = p * L;
    // }
    /* Number of samples for weighted/pool sampling */
    K = p * L;
    Tdiff_pool.resize(K, p);
    if (ctx.get_verbose())
      std::cout << "Size of the genotype pool: " << K << std::endl;
  }

  if (ctx.get_verbose()) {
    std::cout << "Initial value of the error rate - epsilon: "
              << model.get_epsilon() << std::endl;
    std::cout << "Initial value of the rate parameters - lambda: "
              << model.get_lambda().transpose() << std::endl;
  }

  for (unsigned int iter = 0; iter < control_EM.max_iter; ++iter) {

    MatrixXd T_pool;
    if (iter == update_step_size) {
      avg_lambda_current /= control_EM.update_step_size;
      avg_eps_current /= control_EM.update_step_size;
      avg_llhood /= control_EM.update_step_size;
      if (tol_comparison) {
        if (std::abs(avg_eps - avg_eps_current) <= control_EM.tol &&
            ((avg_lambda - avg_lambda_current).array().abs() <= control_EM.tol).all())
          break;
        // L *= 2;
      }
      avg_lambda = avg_lambda_current;
      avg_eps = avg_eps_current;

      update_step_size += control_EM.update_step_size;
      tol_comparison = true;

      /* Restart averaging */
      avg_lambda_current = VectorXd::Zero(p);
      avg_eps_current = 0.0;
      avg_llhood = 0.0;
    }

    /* E step
     * Conditional expectation for the sufficient statistics per observation
     * and event
     */
    if (sampling == "add-remove") {
      scale_cumulative = scale_path_to_mutation(model);
    } else if (sampling == "pool") {
      /* All threads share the same pool of mutation times */
      T_pool.resize(K, p);
      T_pool = sample_times(K, model, Tdiff_pool, ctx.rng);
    }

    N_eff = 0;
    obs_llhood = 0.0;
    expected_dist = 0.0;

    #ifdef _OPENMP
      omp_set_num_threads(thrds);
    #endif
    auto rngs = ctx.get_auxiliary_rngs(thrds);

    #pragma omp parallel for reduction(+:obs_llhood) reduction(+:expected_dist) reduction(+:N_eff) schedule(static)
    for (unsigned int i = 0; i < N; ++i) {
      VectorXi d_pool;
      if (sampling == "pool") {
        VectorXd T_sampling(K);
        if (sampling_times_available)
          T_sampling.setConstant(times(i));

        MatrixXb genotype_pool =
          generate_genotypes(T_pool, model, T_sampling, (*rngs)[omp_get_thread_num()],
                             sampling_times_available);
        d_pool = hamming_dist_mat(genotype_pool, obs.row(i));
      }
      DataImportanceSampling importance_sampling = importance_weight(
        obs.row(i), L, model, times(i), sampling, scale_cumulative, d_pool,
        Tdiff_pool, control_EM.neighborhood_dist, (*rngs)[omp_get_thread_num()],
        sampling_times_available);

      double aux = importance_sampling.w.sum();
      if (aux > 0) {
        /* Only consider observations with at least one feasible sample */
        N_eff += weights(i);
        int L_eff = L;
        if (sampling == "backward")
          L_eff = (importance_sampling.w.array() > 0).count();
        obs_llhood += std::log(weights(i) * aux / L_eff);
        expected_dist += weights(i) *
          importance_sampling.w.dot(importance_sampling.dist.cast<double>()) / aux;
        expected_Tdiff.row(i) =
          (importance_sampling.Tdiff.transpose() * importance_sampling.w) / aux;
      }
    }

    /* M-step */
    model.set_epsilon(expected_dist / (N_eff * p));
    Tdiff_colsum = weights * expected_Tdiff;
    model.set_lambda((Tdiff_colsum / N_eff).array().inverse(), control_EM.max_lambda);

    avg_lambda_current +=  model.get_lambda();
    avg_eps_current += model.get_epsilon();
    avg_llhood += obs_llhood;

    if (iter + 1 == control_EM.max_iter) {
      unsigned int num_iter = control_EM.max_iter - update_step_size +
        control_EM.update_step_size;
      avg_lambda_current /= num_iter;
      avg_eps_current /= num_iter;
      avg_llhood /= num_iter;
    }

    if (ctx.get_verbose()) {
      if (iter == 0)
        std::cout << "llhood\tepsilon\tlambdas" << std::endl;
      std::cout << obs_llhood << "\t" << model.get_epsilon() << "\t"
                << model.get_lambda().transpose() << std::endl;
    }
  }

  model.set_lambda(avg_lambda_current);
  model.set_epsilon(avg_eps_current);
  model.set_llhood(avg_llhood);

  return avg_llhood;
}

RcppExport SEXP _complete_log_likelihood(
    SEXP lambdaSEXP, SEXP epsSEXP, SEXP TdiffSEXP, SEXP distSEXP, SEXP WSEXP) {

  using namespace Rcpp;
  try {
    // Convert input to C++ types
    const MapVecd lambda(as<MapVecd>(lambdaSEXP));
    double eps = as<double>(epsSEXP);
    const MapMatd Tdiff(as<MapMatd>(TdiffSEXP));
    const MapVecd dist(as<MapVecd>(distSEXP));
    float W = as<float>(WSEXP);

    // Call the underlying C++ function
    double res = complete_log_likelihood(lambda, eps, Tdiff, dist, W, false);

    // Return the result as a SEXP
    return wrap( res );
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}

RcppExport SEXP _obs_log_likelihood(
    SEXP obsSEXP, SEXP posetSEXP, SEXP lambdaSEXP, SEXP epsSEXP,
    SEXP weightsSEXP, SEXP timesSEXP, SEXP LSEXP, SEXP samplingSEXP,
    SEXP neighborhood_distSEXP, SEXP lambda_sSEXP,
    SEXP sampling_times_availableSEXP, SEXP thrdsSEXP, SEXP seedSEXP) {

  using namespace Rcpp;
  try {
    // Convert input to C++ types
    const MatrixXb& obs = as<MatrixXb>(obsSEXP);
    const MapMati poset(as<MapMati>(posetSEXP));
    const MapVecd lambda(as<MapVecd>(lambdaSEXP));
    const double eps = as<double>(epsSEXP);
    const MapRowVecd weights(as<MapRowVecd>(weightsSEXP));
    const MapVecd times(as<MapVecd>(timesSEXP));
    const unsigned int L = as<unsigned int>(LSEXP);
    const std::string& sampling = as<std::string>(samplingSEXP);
    const unsigned int neighborhood_dist = as<unsigned int>(neighborhood_distSEXP);
    const float lambda_s = as<float>(lambda_sSEXP);
    const bool sampling_times_available = as<bool>(sampling_times_availableSEXP);
    const int thrds = as<int>(thrdsSEXP);
    const int seed = as<int>(seedSEXP);

    // Call the underlying C++ function
    Context ctx(seed);
    double llhood = obs_log_likelihood(
      obs, poset, lambda, eps, weights, times, L, sampling, neighborhood_dist,
      ctx, lambda_s, sampling_times_available, thrds);

    // Return the result as a SEXP
    return wrap( llhood );
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}

//' @noRd
//' @param update_step_sizeSEXP Evaluate convergence of parameter every
//' 'update_step_size' and increase number of samples, 'L', in order to reach a
//' desirable 'tol'
//' @param tolSEXP Convergence tolerance for rate parameters
RcppExport SEXP _MCEM_hcbn(
    SEXP ilambdaSEXP, SEXP posetSEXP, SEXP obsSEXP, SEXP timesSEXP,
    SEXP lambda_sSEXP, SEXP epsSEXP, SEXP weightsSEXP, SEXP LSEXP,
    SEXP samplingSEXP, SEXP max_iterSEXP, SEXP update_step_sizeSEXP,
    SEXP tolSEXP, SEXP max_lambdaSEXP, SEXP neighborhood_distSEXP,
    SEXP sampling_times_availableSEXP, SEXP thrdsSEXP, SEXP verboseSEXP,
    SEXP seedSEXP) {

  using namespace Rcpp;
  try {
    /* Convert input to C++ types */
    VectorXd ilambda = as<MapVecd>(ilambdaSEXP);
    const MapMati poset(as<MapMati>(posetSEXP));
    const MatrixXb& obs = as<MatrixXb>(obsSEXP);
    const MapVecd times(as<MapVecd>(timesSEXP));
    const float lambda_s = as<float>(lambda_sSEXP);
    const double eps = as<double>(epsSEXP);
    const MapRowVecd weights(as<MapRowVecd>(weightsSEXP));
    unsigned int L = as<unsigned int>(LSEXP);
    const std::string& sampling = as<std::string>(samplingSEXP);
    const unsigned int max_iter = as<unsigned int>(max_iterSEXP);
    const unsigned int update_step_size = as<unsigned int>(update_step_sizeSEXP);
    const double tol = as<double>(tolSEXP);
    const float max_lambda = as<float>(max_lambdaSEXP);
    const unsigned int neighborhood_dist = as<unsigned int>(neighborhood_distSEXP);
    const bool sampling_times_available = as<bool>(sampling_times_availableSEXP);
    const int thrds = as<int>(thrdsSEXP);
    const bool verbose = as<bool>(verboseSEXP);
    const int seed = as<int>(seedSEXP);

    const auto p = poset.rows(); // Number of mutations / events
    edge_container edge_list = adjacency_mat2list(poset);
    Model M(edge_list, p, lambda_s);
    M.set_lambda(ilambda, max_lambda);
    M.set_epsilon(eps);
    M.has_cycles();
    if (M.cycle)
      throw not_acyclic_exception();
    M.topological_sort();

    ControlEM control_EM(max_iter, update_step_size, tol, max_lambda, 
                         neighborhood_dist);

    /* Call the underlying C++ function */
    Context ctx(seed, verbose);
    double llhood = MCEM_hcbn(
      M, obs, times, weights, L, sampling, control_EM,
      sampling_times_available, thrds, ctx);

    /* Return the result as a SEXP */
    return List::create(_["lambda"]=M.get_lambda(), _["eps"]=M.get_epsilon(),
                        _["llhood"]=llhood);
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}

RcppExport SEXP _importance_weight_genotype(
    SEXP genotypeSEXP, SEXP LSEXP, SEXP posetSEXP, SEXP lambdaSEXP,
    SEXP epsSEXP, SEXP timeSEXP, SEXP samplingSEXP, SEXP scale_cumulativeSEXP,
    SEXP d_poolSEXP, SEXP Tdiff_poolSEXP, SEXP neighborhood_distSEXP,
    SEXP lambda_sSEXP, SEXP sampling_times_availableSEXP, SEXP seedSEXP) {

  using namespace Rcpp;
  try {
    /* Convert input to C++ types */
    const RowVectorXb& genotype = as<RowVectorXb>(genotypeSEXP);
    const unsigned int L = as<unsigned int>(LSEXP);
    const MapMati poset(as<MapMati>(posetSEXP));
    const MapVecd lambda(as<MapVecd>(lambdaSEXP));
    const double eps = as<double>(epsSEXP);
    const double time = as<double>(timeSEXP);
    const std::string& sampling = as<std::string>(samplingSEXP);
    const MapVecd scale_cumulative(as<MapVecd>(scale_cumulativeSEXP));
    const MapVeci d_pool(as<MapVeci>(d_poolSEXP));
    const MapMatd Tdiff_pool(as<MapMatd>(Tdiff_poolSEXP));
    const unsigned int neighborhood_dist = as<unsigned int>(neighborhood_distSEXP);
    const float lambda_s = as<float>(lambda_sSEXP);
    const bool sampling_times_available = as<bool>(sampling_times_availableSEXP);
    const int seed = as<int>(seedSEXP);

    const auto p = poset.rows(); // Number of mutations / events
    edge_container edge_list = adjacency_mat2list(poset);
    Model M(edge_list, p, lambda_s);
    M.set_lambda(lambda);
    M.set_epsilon(eps);
    M.has_cycles();
    if (M.cycle)
      throw not_acyclic_exception();
    M.topological_sort();

    /* Call the underlying C++ function */
    Context ctx(seed);
    DataImportanceSampling w = importance_weight(
      genotype, L, M, time, sampling, scale_cumulative, d_pool, Tdiff_pool,
      neighborhood_dist, ctx.rng, sampling_times_available);

    /* Return the result as a SEXP */
    return List::create(_["w"]=w.w, _["dist"]=w.dist, _["Tdiff"]=w.Tdiff);
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}

RcppExport SEXP _importance_weight(
    SEXP obsSEXP, SEXP LSEXP, SEXP posetSEXP, SEXP lambdaSEXP,
    SEXP epsSEXP, SEXP timesSEXP, SEXP samplingSEXP, SEXP neighborhood_distSEXP,
    SEXP lambda_sSEXP, SEXP sampling_times_availableSEXP, SEXP thrdsSEXP,
    SEXP seedSEXP) {

  using namespace Rcpp;
  try {
    /* Convert input to C++ types */
    const MatrixXb& obs = as<MatrixXb>(obsSEXP);
    const unsigned int L = as<unsigned int>(LSEXP);
    const MapMati poset(as<MapMati>(posetSEXP));
    const MapVecd lambda(as<MapVecd>(lambdaSEXP));
    const double eps = as<double>(epsSEXP);
    const MapVecd times(as<MapVecd>(timesSEXP));
    const std::string& sampling = as<std::string>(samplingSEXP);
    const unsigned int neighborhood_dist = as<unsigned int>(neighborhood_distSEXP);
    const float lambda_s = as<float>(lambda_sSEXP);
    const bool sampling_times_available = as<bool>(sampling_times_availableSEXP);
    const int thrds = as<int>(thrdsSEXP);
    const int seed = as<int>(seedSEXP);

    const unsigned int N = obs.rows(); // Number of observations / genotypes
    const auto p = poset.rows();       // Number of mutations / events
    VectorXd scale_cumulative;
    VectorXd w_sum(N);
    VectorXd w_sum_sqrt(N);
    VectorXd expected_dist(N);
    MatrixXd expected_Tdiff(N, p);
    edge_container edge_list = adjacency_mat2list(poset);
    Model M(edge_list, p, lambda_s);
    M.set_lambda(lambda);
    M.set_epsilon(eps);
    M.has_cycles();
    if (M.cycle)
      throw not_acyclic_exception();
    M.topological_sort();

    Context ctx(seed);
    unsigned int K = 0;
    MatrixXd Tdiff_pool;
    MatrixXd T_pool;
    std::vector<int> L_eff(N);
    if (sampling == "add-remove") {
      scale_cumulative.resize(p);
      scale_cumulative = scale_path_to_mutation(M);
    } else if (sampling == "pool") {
      K = p * L;
      Tdiff_pool.resize(K, p);
      T_pool.resize(K, p);
      T_pool = sample_times(K, M, Tdiff_pool, ctx.rng);
    }

    #ifdef _OPENMP
    omp_set_num_threads(thrds);
    #endif
    auto rngs = ctx.get_auxiliary_rngs(thrds);

    #pragma omp parallel for schedule(static)
    for (unsigned int i = 0; i < N; ++i) {
      VectorXi d_pool;
      if (sampling == "pool") {
        VectorXd T_sampling(K);
        if (sampling_times_available)
          T_sampling.setConstant(times[i]);

        MatrixXb genotype_pool =
          generate_genotypes(T_pool, M, T_sampling, (*rngs)[omp_get_thread_num()],
                             sampling_times_available);
        d_pool = hamming_dist_mat(genotype_pool, obs.row(i));
      }
      /* Call the underlying C++ function */
      DataImportanceSampling w = importance_weight(
        obs.row(i), L, M, times(i), sampling, scale_cumulative, d_pool,
        Tdiff_pool, neighborhood_dist, (*rngs)[omp_get_thread_num()],
        sampling_times_available);

      if (sampling == "backward" || sampling == "bernoulli")
        L_eff[i] = (w.w.array() > 0).count();
      w_sum[i] = w.w.sum();
      w_sum_sqrt[i] = w.w.dot(w.w);
      expected_dist[i] = w.w.dot(w.dist.cast<double>()) / w_sum[i];
      expected_Tdiff.row(i) = (w.Tdiff.transpose() * w.w) / w_sum[i];
    }

    /* Return the result as a SEXP */
    if (sampling == "backward" || sampling == "bernoulli")
      return List::create(_["w"]=w_sum, _["w_sqrt"]=w_sum_sqrt, _["L_eff"]=L_eff,
                          _["dist"]=expected_dist, _["Tdiff"]=expected_Tdiff);
    else
      return List::create(_["w"]=w_sum, _["w_sqrt"]=w_sum_sqrt,
                          _["dist"]=expected_dist, _["Tdiff"]=expected_Tdiff);
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}

RcppExport SEXP _sample_genotypes(
    SEXP NSEXP, SEXP posetSEXP, SEXP lambdaSEXP, SEXP T_eventsSEXP,
    SEXP T_samplingSEXP, SEXP lambda_sSEXP, SEXP sampling_times_availableSEXP,
    SEXP seedSEXP) {

  using namespace Rcpp;
  try {
    /* Convert input to C++ types */
    const unsigned int N = as<unsigned int>(NSEXP);
    const MapMati poset(as<MapMati>(posetSEXP));
    const MapVecd lambda(as<MapVecd>(lambdaSEXP));
    MatrixXd T_events = as<MapMatd>(T_eventsSEXP);
    VectorXd T_sampling = as<MapVecd>(T_samplingSEXP);
    const float lambda_s = as<float>(lambda_sSEXP);
    const bool sampling_times_available = as<bool>(sampling_times_availableSEXP);
    const int seed = as<int>(seedSEXP);

    const auto p = poset.rows(); // Number of mutations / events
    edge_container edge_list = adjacency_mat2list(poset);
    Model M(edge_list, p, lambda_s);
    M.set_lambda(lambda);
    M.has_cycles();
    if (M.cycle)
      throw not_acyclic_exception();
    M.topological_sort();

    /* Call the underlying C++ function */
    Context ctx(seed);
    MatrixXb samples = sample_genotypes(N, M, T_events, T_sampling, ctx.rng,
                                        sampling_times_available);

    /* Return the result as a SEXP*/
    return List::create(_["samples"]=samples, _["Tdiff"]=T_events,
                        _["sampling_time"]=T_sampling);
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}

RcppExport SEXP _sample_times(
    SEXP NSEXP, SEXP posetSEXP, SEXP lambdaSEXP, SEXP T_eventsSEXP,
    SEXP seedSEXP) {

  using namespace Rcpp;
  try {
    /* Convert input to C++ types */
    const unsigned int N = as<unsigned int>(NSEXP);
    const MapMati poset(as<MapMati>(posetSEXP));
    const MapVecd lambda(as<MapVecd>(lambdaSEXP));
    MatrixXd T_events = as<MapMatd>(T_eventsSEXP);
    const int seed = as<int>(seedSEXP);

    const auto p = poset.rows(); // Number of mutations / events
    edge_container edge_list = adjacency_mat2list(poset);
    Model M(edge_list, p);
    M.set_lambda(lambda);
    M.has_cycles();
    if (M.cycle)
      throw not_acyclic_exception();
    M.topological_sort();

    /* Call the underlying C++ function */
    Context ctx(seed);
    MatrixXd times = sample_times(N, M, T_events, ctx.rng);

    /* Return the result as a SEXP */
    return List::create(_["mutation_times"]=times, _["Tdiff"]=T_events);
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}

RcppExport SEXP _generate_genotypes(
    SEXP T_events_sumSEXP, SEXP posetSEXP, SEXP T_samplingSEXP,
    SEXP lambda_sSEXP, SEXP sampling_times_availableSEXP, SEXP seedSEXP) {

  using namespace Rcpp;
  try {
    /* Convert input to C++ types */
    const MapMatd T_events_sum(as<MapMatd>(T_events_sumSEXP));
    const MapMati poset(as<MapMati>(posetSEXP));
    VectorXd T_sampling = as<MapVecd>(T_samplingSEXP);
    const float lambda_s = as<float>(lambda_sSEXP);
    const bool sampling_times_available = as<bool>(sampling_times_availableSEXP);
    const int seed = as<int>(seedSEXP);

    const auto p = poset.rows(); // Number of mutations / events
    edge_container edge_list = adjacency_mat2list(poset);
    Model M(edge_list, p, lambda_s);
    M.has_cycles();
    if (M.cycle)
      throw not_acyclic_exception();
    M.topological_sort();

    Context ctx(seed);

    /* Call the underlying C++ function */
    MatrixXb samples = generate_genotypes(T_events_sum, M, T_sampling, ctx.rng,
                                          sampling_times_available);

    /* Return the result as a SEXP */
    return List::create(_["samples"]=samples, _["sampling_time"]=T_sampling);
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}
