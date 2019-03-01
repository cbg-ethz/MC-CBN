/** mccbn: large-scale inference on conjunctive Bayesian networks
 *  MCEM for the hidden conjuctive Bayesian network model
 *
 * @author Susana Posada Céspedes
 * @email susana.posada@bsse.ethz.ch
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "mcem.h"
#include <random>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
typedef Map<MatrixXd> MapMatd;

class Initializer {
public:
  Initializer() {
    Eigen::initParallel();
  }
} initializer;

// void set_seed(std::mt19937& rng, int seed) {
//   std::ranlux24_base seeder_rng = std::ranlux24_base(seed);
//
//   // adapted from https://stackoverflow.com/a/15509942
//   std::array<int, std::mt19937::state_size> seed_data;
//   std::generate_n(seed_data.data(), seed_data.size(), std::ref(seeder_rng));
//   std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
//
//   rng.seed(seq);
// }

// void seed_rng(std::mt19937& rng, int seed) {
//   set_seed(rng, seed);
// }


//' @param N number of samples to be drawn
//' @param lambda rate
VectorXd rexp_std(const unsigned int N, const double lambda,
                  Context::rng_type& rng) {
  std::exponential_distribution<double> distribution(lambda);
  VectorXd T(N);
  
  for (unsigned int i = 0; i < N; ++i)
    T[i] = distribution(rng);
  
  return T;
}

std::vector<int> rdiscrete_std(const unsigned int N, const VectorXd weights,
                               std::mt19937& rng) {
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
        log_prob[i] = log(eps + DBL_EPSILON) * dist[i] +
          log(1 - eps - DBL_EPSILON) * (p - dist[i]);
      } else {
        log_prob[i] = 0;
      }
    }
  } else {
    log_prob = (log(eps) * dist).array() + log(1 - eps) * (p - dist.array());
  }
  return log_prob;
}

//' Compute complete-data log-likelihood or (equivalently) hidden log-likelihood
double complete_log_likelihood(const VectorXd& lambda, const double eps,
                               const MatrixXd& Tdiff, const VectorXd& dist,
                               const float W) {

  const unsigned int p = lambda.size();
  const unsigned int N = dist.size();
  double llhood;

  llhood = W * lambda.array().log().sum() - (Tdiff * lambda).sum();
  if (eps == 0) {
    for (unsigned int i = 0; i < N; ++i) {
      if (dist(i) != 0) {
        llhood += log(eps + DBL_EPSILON) * dist(i) +
          log(1 - eps - DBL_EPSILON) * (p - dist(i));
      }
    }
  } else {
    llhood += log(eps) * dist.sum() + log(1 - eps) * (p - dist.array()).sum();
  }
  return llhood;
}

//' Compute observed log-likelihood
double obs_log_likelihood(
    const MatrixXb& obs, const MatrixXi& poset, const VectorXd& lambda,
    const VectorXi& topo_path, const double eps, const VectorXd& times,
    const unsigned int L, const std::string& sampling, const unsigned int version,
    const float perturb_prob, const MatrixXb& genotype_pool,
    const MatrixXd& Tdiff_pool, Context& ctx, const float lambda_s=1.0,
    const bool sampling_times_available=false, const unsigned int thrds=1) {

  const unsigned int p = poset.rows(); // Number of mutations / events
  const unsigned int N = obs.rows();   // Number of observations / genotypes
  double llhood = 0;
  
  Model model(p, lambda_s);
  model.fill_poset(poset);
  model.parents();
  model.set_lambda(lambda);
  model.set_epsilon(eps);

  #ifdef _OPENMP
    omp_set_num_threads(thrds);
  #endif
  auto rngs = ctx.get_auxiliary_rngs(thrds);

  #pragma omp parallel for reduction(+:llhood) schedule(static)
  for (unsigned int i = 0; i < N; ++i) {
    VectorXi dist_pool;
    if (sampling == "backward")
      dist_pool = hamming_dist_mat(genotype_pool, obs.row(i));

    DataImportanceSampling importance_sampling = importance_weight(
      obs.row(i), L, model, topo_path, times[i], sampling, version,
      perturb_prob, dist_pool, Tdiff_pool, rngs[omp_get_thread_num()],
      sampling_times_available);
    llhood = llhood + std::log(importance_sampling.w.sum() / L);
  }
  return llhood;
}

MatrixXb sample_genotypes(
    const unsigned int N, const Model& model, const VectorXi& topo_path,
    MatrixXd& T_events, VectorXd& T_sampling, Context::rng_type& rng,
    const bool sampling_times_available=false) {
  
  // Initialization and instantiation of variables
  const unsigned int p = model.size();     // Number of mutations / events
  MatrixXd T_sum_events(N, p);
  MatrixXb obs;
  obs.setZero(N, p);
  
  // Generate occurence times T_events_{j} ~ Exp(lambda_{j})
  for (unsigned int j = 0; j < p; ++j)
    T_events.col(j) = rexp_std(N, model.get_lambda(j), rng);
    // NumericVector aux = rexp(N, lambda[j]);
    // T_events.col(j) = as<VectorXd>(aux);
  
  // Use sampling times when available
  if (!sampling_times_available)
    T_sampling = rexp_std(N, model.get_lambda_s(), rng);
    // NumericVector aux = rexp(N, model.get_lambda_s());
    // T_sampling = as<VectorXd>(aux);
  
  // Rcout << T_sampling << std::endl; // DBG T_sampling
  
  for (unsigned int k = 0; k < p; ++k) {
    int j = topo_path(k);
    for (unsigned int i = 0; i < N; ++i) {
      double T_max = 0.0;
      for (int u = 0; u < model.N_pa[j]; ++u)
        if (T_sum_events(i, model.pa[j][u]) > T_max)
          T_max = T_sum_events(i, model.pa[j][u]);
      T_sum_events(i, j) = T_events(i, j) + T_max;
      if (T_sum_events(i, j) <= T_sampling(i))
        obs(i, j) = 1; 
    }
  }
  
  return obs;
}

int hamming_dist(const VectorXi& x, const VectorXi& y) {
  return (x - y).array().abs().sum();
}

VectorXi hamming_dist_mat(const MatrixXb& x, const RowVectorXb& y) {
  const int N = x.rows();
  return (x.array() != y.replicate(N, 1).array()).rowwise().count().cast<int>();
  // return (x.rowwise() - y).array().abs().rowwise().sum();
}

//' Compute importance weights and sufficient statistics by sampling
//'
//' @noRd
//' @param genotype
//' @param L number of samples
//' @param model
//' @param topo_path
//' @param time
//' @param K number of samples for weighted sampling
//' @return  returns importance weights and sufficient statistics
DataImportanceSampling importance_weight(
    const RowVectorXb& genotype, const unsigned int L, const Model& model,
    const VectorXi& topo_path, const double time, const std::string& sampling,
    const unsigned int version, const float perturb_prob,
    const VectorXi& dist_pool, const MatrixXd& Tdiff_pool,
    Context::rng_type& rng, const bool sampling_times_available=false) {
  
  // Initialization and instantiation of variables
  const unsigned int p = model.size(); // Number of mutations / events
  MatrixXb samples(L, p);
  DataImportanceSampling importance_sampling(L, p);
  
  if (sampling == "naive") {
    /* Generate L samples from poset with parameters lambda and lambda_s.
     * In particular, epsilon is zero (default value) - because the idea is to
     * generate samples of X (true genotype)
     */
    VectorXd T_sampling(L);
    if (sampling_times_available)
      T_sampling.setConstant(time);

    samples = sample_genotypes(
      L, model, topo_path, importance_sampling.Tdiff, T_sampling, rng,
      sampling_times_available);
    importance_sampling.dist = hamming_dist_mat(samples, genotype);
    VectorXd d = importance_sampling.dist.cast<double>();
    importance_sampling.w = pow(model.get_epsilon(), d.array()) *
      pow(1 - model.get_epsilon(), p - d.array());
  } else if (sampling == "add-remove") {
    //TODO
  } else if (sampling == "backward") {
    const unsigned int K = dist_pool.size();
    VectorXd log_prob_Y_X(L);
    VectorXd log_proposal(L);
    VectorXd q_prob_selected(L);
    
    VectorXd d_pool = dist_pool.cast<double>();
    VectorXd q_prob = pow(model.get_epsilon(), d_pool.array()) *
      pow(1 - model.get_epsilon(), p - d_pool.array());
    q_prob /= q_prob.sum();
    
    // Draw L samples with replacement and with weights q_prob
    std::vector<int> idxs_sample = rdiscrete_std(L, q_prob, rng);
    int idx; 
    for (unsigned int l = 0; l < L; ++l) {
      // TODO: can this be done better?
      idx = idxs_sample[l];
      importance_sampling.dist(l) = dist_pool(idx);
      q_prob_selected(l) = q_prob(idx);
      importance_sampling.Tdiff.row(l) = Tdiff_pool.row(idx);
    }

    log_prob_Y_X = log_bernoulli_process(importance_sampling.dist.cast<double>(),
                                         model.get_epsilon(), p);
    
    // Computing density of the proposal - for correction
    log_proposal = q_prob_selected.array().log() + log(K);
    
    importance_sampling.w = (log_prob_Y_X - log_proposal).array().exp();
  }

  return importance_sampling;
}

//' Compute importance weights and sufficient statistics by sampling
//'
//' @noRd
double MCEM_hcbn(
    Model& model, const MatrixXb& obs, const VectorXd& times,
    const VectorXi& topo_path, const RowVectorXd& weights,
    const unsigned int L, const std::string& sampling,
    const unsigned int version, const float perturb_prob,
    const unsigned int max_iter, const float burn_in, const float max_lambda,
    const bool sampling_times_available, const unsigned int thrds,
    Context& ctx) {
  
  // Initialization and instantiation of variables
  const unsigned int p = obs.cols();  // Number of mutations / events
  const unsigned int N = obs.rows();  // Number of observations / genotypes
  float W = weights.sum();            // Number of (weighted) observations
  unsigned int K = 0;
  VectorXd avg_lambda = VectorXd::Zero(p);
  double avg_eps = 0, avg_llhood = 0, llhood = 0;
  VectorXd expected_dist(N);
  MatrixXd expected_Tdiff(N, p);
  VectorXd Tdiff_colsum(p);
  
  // Model model(p);
  // model.fill_poset(poset);
  // model.parents();
  
  // MatrixXb genotype_pool;
  // MatrixXd Tdiff_pool;
  if (sampling == "backward") {
    // if (p < 19) {
    //   K = std::max((unsigned int) pow(2, p + 1), 2 * L);
    // } else {
    //   // K = 100000000 / (sizeof(double) * p); //NOTE:empirically 100000 limit on runtime
    //   K = p * L;
    // }
    K = p * L;
    // genotype_pool.resize(K, p);
    // Tdiff_pool.resize(K, p);
    std::cout << "Size of the genotype pool: " << K << std::endl;
  }

  if (ctx.get_verbose())
    std::cout << "Initial value for rate parameters - lambda: "
              << model.get_lambda().transpose() << std::endl;
  
  const unsigned int record_iter = std::max(int(burn_in * max_iter), 1); 
  
  for (unsigned int iter = 0; iter < max_iter; ++iter) {
    
    /* E step
     * Conditional expectation for the sufficient statistics per observation
     * and event
     */
    // NOTE: All threads sharing same pool
    // if (sampling == "backward") {
    //   VectorXd T_sampling;
    //   if (sampling_times_available) {
    //     T_sampling.resize(K);
    //     NumericVector times_aux(wrap(times));
    //     NumericVector aux = sample(times_aux, K, true);
    //     T_sampling = as<MapVecd>(aux);
    //   }
    //   genotype_pool = sample_genotypes(
    //     K, model, topo_path, Tdiff_pool, T_sampling, ctx.rng,
    //     sampling_times_available);
    // }
    #ifdef _OPENMP
      omp_set_num_threads(thrds);
    #endif
    auto rngs = ctx.get_auxiliary_rngs(thrds);
    // std::vector<std::mt19937> rngs(thrds);
    // for (unsigned int t = 0; t < thrds; ++t)
    //   set_seed(rngs[t], ctx.rng());
      // set_seed(rngs[t], main_rng());
    #pragma omp parallel for schedule(static)
    for (unsigned int i = 0; i < N; ++i) {
      VectorXi d_pool;
      MatrixXd Tdiff_pool;
      if (sampling == "backward") {
        Tdiff_pool.resize(K, p);
        VectorXd T_sampling;
        if (sampling_times_available) {
          T_sampling.resize(K);
          T_sampling.setConstant(times(i));
          // Sample with replacement from vector containing sampling times
          // TODO
          // NumericVector times_aux(wrap(times));
          // NumericVector aux = sample(times_aux, K, true);
          // T_sampling = as<MapVecd>(aux);
        }
        MatrixXb genotype_pool = sample_genotypes(
          K, model, topo_path, Tdiff_pool, T_sampling,
          rngs[omp_get_thread_num()], sampling_times_available);
        d_pool = hamming_dist_mat(genotype_pool, obs.row(i));
      }
      DataImportanceSampling importance_sampling = importance_weight(
        obs.row(i), L, model, topo_path, times(i), sampling, version,
        perturb_prob, d_pool, Tdiff_pool, rngs[omp_get_thread_num()],
        sampling_times_available);

      expected_dist(i) =
        importance_sampling.w.dot(importance_sampling.dist.cast<double>()) /
          importance_sampling.w.sum();
      expected_Tdiff.row(i) = 
        (importance_sampling.Tdiff.transpose() * importance_sampling.w) /
        importance_sampling.w.sum();
    }
    
    // M-step
    model.set_epsilon(expected_dist.sum() / (N * p));
    Tdiff_colsum = weights * expected_Tdiff;
    model.set_lambda((Tdiff_colsum / W).array().inverse(), max_lambda);
    
    llhood = complete_log_likelihood(
      model.get_lambda(), model.get_epsilon(), expected_Tdiff, expected_dist,
      W);
    if (ctx.get_verbose()) {
      if (iter == 0)
        std::cout << "llhood\tepsilon\tlambdas" << std::endl;
      std::cout << llhood << "\t" << model.get_epsilon() << "\t"
                << model.get_lambda().transpose() << std::endl;
    }

    if (iter >= record_iter) {
      avg_lambda +=  model.get_lambda();
      avg_eps += model.get_epsilon();
      avg_llhood += llhood;
    }
  }
  
  avg_lambda /= (max_iter - record_iter);
  avg_eps /= (max_iter - record_iter);
  avg_llhood /= (max_iter - record_iter);
  
  model.set_lambda(avg_lambda);
  model.set_epsilon(avg_eps);
  model.set_llhood(avg_llhood);
  
  return avg_llhood;
}

RcppExport SEXP _complete_log_likelihood(
    SEXP lambdaSEXP, SEXP epsSEXP, SEXP TdiffSEXP, SEXP distSEXP, SEXP WSEXP) {
  
  // Convert input to C++ types
  const MapVecd lambda(as<MapVecd>(lambdaSEXP));
  double eps = as<double>(epsSEXP);
  const MapMatd Tdiff(as<MapMatd>(TdiffSEXP));
  const MapVecd dist(as<MapVecd>(distSEXP));
  float W = as<float>(WSEXP);

  // Call the underlying C++ function
  double res = complete_log_likelihood(lambda, eps, Tdiff, dist, W);

  // Return the result as a SEXP
  return wrap( res );
}

RcppExport SEXP _obs_log_likelihood(
    SEXP obsSEXP, SEXP posetSEXP, SEXP lambdaSEXP, SEXP topo_pathSEXP,
    SEXP epsSEXP, SEXP timesSEXP, SEXP LSEXP, SEXP samplingSEXP, SEXP versionSEXP,
    SEXP perturb_probSEXP, SEXP genotype_poolSEXP, //SEXP dist_poolSEXP,
    SEXP Tdiff_poolSEXP, SEXP lambda_sSEXP,  SEXP sampling_times_availableSEXP,
    SEXP thrdsSEXP, SEXP seedSEXP) {
  
  try {
    // Convert input to C++ types
    const MatrixXb& obs = as<MatrixXb>(obsSEXP);
    const MapMati poset(as<MapMati>(posetSEXP));
    const MapVecd lambda(as<MapVecd>(lambdaSEXP));
    const MapVeci topo_path(as<MapVeci>(topo_pathSEXP));
    const double eps = as<double>(epsSEXP);
    const MapVecd times(as<MapVecd>(timesSEXP));
    const unsigned int L = as<unsigned int>(LSEXP);
    const std::string& sampling = as<std::string>(samplingSEXP);
    const unsigned int version = as<unsigned int>(versionSEXP);
    const float perturb_prob = as<float>(perturb_probSEXP);
    // const MapMati dist_pool(as<MapMati>(dist_poolSEXP));
    const MatrixXb& genotype_pool = as<MatrixXb>(genotype_poolSEXP);
    const MapMatd Tdiff_pool(as<MapMatd>(Tdiff_poolSEXP));
    const float lambda_s = as<float>(lambda_sSEXP);
    const bool sampling_times_available = as<bool>(sampling_times_availableSEXP);
    const int thrds = as<int>(thrdsSEXP);
    const int seed = as<int>(seedSEXP);
    
    // Call the underlying C++ function
    // seed_rng(seed);
    Context ctx(seed);
    double llhood = obs_log_likelihood(
      obs, poset, lambda, topo_path, eps, times, L, sampling, version,
      perturb_prob, genotype_pool, Tdiff_pool, ctx, lambda_s,
      sampling_times_available, thrds);
    
    // Return the result as a SEXP
    return wrap( llhood );
  } catch (std::exception &ex) {
    // NOTE: reference to 'exception' is ambiguous, 'std::' required
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue;
}

RcppExport SEXP _MCEM_hcbn(
    SEXP ilambdaSEXP, SEXP posetSEXP, SEXP obsSEXP, SEXP timesSEXP, 
    SEXP lambda_sSEXP, SEXP topo_pathSEXP, SEXP epsSEXP, SEXP weightsSEXP,
    SEXP LSEXP, SEXP samplingSEXP, SEXP versionSEXP, SEXP perturb_probSEXP,
    SEXP max_iterSEXP, SEXP burn_inSEXP, SEXP max_lambdaSEXP,
    SEXP sampling_times_availableSEXP, SEXP thrdsSEXP, SEXP verboseSEXP,
    SEXP seedSEXP) { 

  // Convert input to C++ types
  VectorXd ilambda = as<MapVecd>(ilambdaSEXP);
  const MapMati poset(as<MapMati>(posetSEXP)); 
  const MatrixXb& obs = as<MatrixXb>(obsSEXP);
  const MapVecd times(as<MapVecd>(timesSEXP));
  const float lambda_s = as<float>(lambda_sSEXP);
  const MapVeci topo_path(as<MapVeci>(topo_pathSEXP));
  const double eps = as<double>(epsSEXP);
  const MapRowVecd weights(as<MapRowVecd>(weightsSEXP));
  const unsigned int L = as<unsigned int>(LSEXP);
  const std::string& sampling = as<std::string>(samplingSEXP);
  const unsigned int version = as<unsigned int>(versionSEXP);
  const float perturb_prob = as<float>(perturb_probSEXP);
  const unsigned int max_iter = as<unsigned int>(max_iterSEXP);
  const float burn_in = as<float>(burn_inSEXP);
  const float max_lambda = as<float>(max_lambdaSEXP);
  const bool sampling_times_available = as<bool>(sampling_times_availableSEXP);
  const int thrds = as<int>(thrdsSEXP);
  const bool verbose = as<bool>(verboseSEXP);
  const int seed = as<int>(seedSEXP);
  
  const unsigned int p = poset.rows();
  Model M(p, lambda_s);
  M.fill_poset(poset);
  M.parents();
  M.set_lambda(ilambda);
  M.set_epsilon(eps);

  // Call the underlying C++ function
  Context ctx(seed, verbose);
  double llhood = MCEM_hcbn(
    M, obs, times, topo_path, weights, L, sampling, version, perturb_prob,
    max_iter, burn_in, max_lambda, sampling_times_available, thrds, ctx);
    
  // Return the result as a SEXP
  return List::create(_["lambda"]=M.get_lambda(), _["eps"]=M.get_epsilon(),
                      _["llhood"]=llhood);
}

RcppExport SEXP _importance_weight(
    SEXP genotypeSEXP, SEXP LSEXP, SEXP posetSEXP, SEXP lambdaSEXP,
    SEXP topo_pathSEXP, SEXP epsSEXP, SEXP timeSEXP, SEXP samplingSEXP, 
    SEXP versionSEXP, SEXP perturb_probSEXP, SEXP d_poolSEXP, 
    SEXP Tdiff_poolSEXP, SEXP lambda_sSEXP, SEXP sampling_times_availableSEXP,
    SEXP seedSEXP) {
  
  try {
    // Convert input to C++ types
    const RowVectorXb& genotype = as<RowVectorXb>(genotypeSEXP);
    const unsigned int L = as<unsigned int>(LSEXP);
    const MapMati poset(as<MapMati>(posetSEXP));
    const MapVecd lambda(as<MapVecd>(lambdaSEXP));
    const MapVeci topo_path(as<MapVeci>(topo_pathSEXP));
    const double eps = as<double>(epsSEXP);
    const double time = as<double>(timeSEXP);
    const std::string& sampling = as<std::string>(samplingSEXP);
    const unsigned int version = as<unsigned int>(versionSEXP);
    const float perturb_prob = as<float>(perturb_probSEXP);
    const MapVeci d_pool(as<MapVeci>(d_poolSEXP)); 
    const MapMatd Tdiff_pool(as<MapMatd>(Tdiff_poolSEXP));
    const float lambda_s = as<float>(lambda_sSEXP);
    const bool sampling_times_available = as<bool>(sampling_times_availableSEXP);
    const int seed = as<int>(seedSEXP);

    const unsigned int p = poset.rows();
    Model M(p, lambda_s);
    M.fill_poset(poset);
    M.parents();
    M.set_lambda(lambda);
    M.set_epsilon(eps);
      
    // Call the underlying C++ function
    Context ctx(seed);
    DataImportanceSampling w = importance_weight(
      genotype, L, M, topo_path, time, sampling, version, perturb_prob, d_pool,
      Tdiff_pool, ctx.rng, sampling_times_available);

    // Return the result as a SEXP
    return List::create(_["w"]=w.w, _["dist"]=w.dist, _["Tdiff"]=w.Tdiff);
  } catch (std::exception &ex) {
    // NOTE: reference to 'exception' is ambiguous, 'std::' required
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue;
}

RcppExport SEXP _sample_genotypes(
    SEXP NSEXP, SEXP posetSEXP, SEXP lambdaSEXP, SEXP topo_pathSEXP,
    SEXP T_eventsSEXP, SEXP T_samplingSEXP, SEXP lambda_sSEXP, 
    SEXP sampling_times_availableSEXP, SEXP seedSEXP) {
  
  // Convert input to C++ types
  const unsigned int N = as<unsigned int>(NSEXP);
  const MapMati poset(as<MapMati>(posetSEXP));
  const MapVecd lambda(as<MapVecd>(lambdaSEXP));
  const MapVeci topo_path(as<MapVeci>(topo_pathSEXP));
  MatrixXd T_events = as<MapMatd>(T_eventsSEXP);
  VectorXd T_sampling = as<MapVecd>(T_samplingSEXP);
  const float lambda_s = as<float>(lambda_sSEXP);
  const bool sampling_times_available = as<bool>(sampling_times_availableSEXP);
  const int seed = as<int>(seedSEXP);
  
  const unsigned int p = poset.rows();
  Model M(p, lambda_s);
  M.fill_poset(poset);
  M.parents();
  M.set_lambda(lambda);
  
  // Rcout << "Before function call " << T_sampling.transpose() << std::endl; // DBG T_sampling
  // Call the underlying C++ function
  Context ctx(seed);
  MatrixXb samples = sample_genotypes(
    N, M, topo_path, T_events, T_sampling, ctx.rng, sampling_times_available);
  // Rcout << "After function call " << T_sampling.transpose() << std::endl; // DBG T_sampling
  
  // Return the result as a SEXP
  return List::create(_["samples"]=samples, _["Tdiff"]=T_events, 
                      _["T_sampling"]=T_sampling);
}