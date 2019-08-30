/** 
 * Network learning using adaptive simulated annealing
 *
 * This file is part of the mccbn package
 * 
 * @author Susana Posada Céspedes
 * @email susana.posada@bsse.ethz.ch
 */

#include <random>
#include <vector>
#include <queue>
#include <iostream>
#include <fstream>
#include <string>

#include <Rcpp.h>
#include <RcppEigen.h>

#include <boost/graph/graph_traits.hpp>
// #include <boost/filesystem.hpp>

#include "mcem.hpp"
#include "not_acyclic_exception.hpp"

using namespace Rcpp;

//' Propose moves for the simulated annealing
//'
//' @noRd
double propose_edge(
    Model& M, const MatrixXb& obs, const double llhood_current,
    double& fraction_compatible, const VectorXd& times,
    const RowVectorXd& weights, const float T, int& num_accept,
    const float factor_fraction_compatible, const unsigned int L,
    const std::string& sampling, const unsigned int version,
    const ControlEM& control_EM, const bool sampling_times_available,
    const unsigned int thrds, Context& ctx) {
  
  const vertices_size_type p = M.size();  // Number of mutations / events
  const auto N = obs.rows();              // Number of observations / genotypes
  
  Model M_new(M);
  double llhood_new = 0.0;
  double fraction_compatible_new = 0.0;
  std::vector<int> idxs_pool(p * p);
  unsigned int v1;
  unsigned int v2;
  
  /* Shuffle pool of indices. Use to draw edges/cover relations */
  std::iota(idxs_pool.begin(), idxs_pool.end(), 0);
  std::shuffle(idxs_pool.begin(), idxs_pool.end(), ctx.rng);
  
  std::uniform_real_distribution<> rand(0.0, 1.0);
  
  /* Propose an update */
  for (unsigned int i = 0; i < p * p; ++i) {
    /* Every time a new move is proposed, start with initial poset */
    M_new = M;
    
    /* Draw one relation at random */
    v1 = (int) idxs_pool[i] / p;
    v2 = idxs_pool[i] % p;
    if (ctx.get_verbose())
      std::cout << "Testing edge: " << v1 << "->" << v2 << std::endl;
    
    if (v1 == v2)
      continue;
    
    /* Change edge
     * if the edge already existed, remove it. Otherwise, add it
     */
    std::pair<boost::graph_traits<Poset>::edge_descriptor, bool> 
      add = boost::add_edge(v1, v2, M_new.poset);
    if (add.second) {
      /* Add edge - Remove the edges that cause the resulting graph not to be
       * transitively reduced
       */
      M_new.has_cycles();
      if (M_new.cycle) {
        if (ctx.get_verbose())
          std::cout << "Cycle!" << std::endl;
        continue;
      } else {
        M_new.topological_sort();
        if (M_new.transitive_reduction_dag()) {
          if (ctx.get_verbose())
            std::cout << "Adding new edge: " <<  add.first << std::endl;
        } else {
          continue;
        }
      }
    } else {
      /* Remove edge*/
      remove_edge(v1, v2, M_new.poset);
      if (ctx.get_verbose())
        std::cout << "Removing edge: " <<  add.first << std::endl;
    }
    
    /* Heuristic based on the fraction of compatible observations with the
     * proposed poset
     */
    fraction_compatible_new =
      (double) num_compatible_observations(obs, M_new, N) / N;
    if (ctx.get_verbose())
      std::cout << "Fraction of compatible observations: "
                << fraction_compatible_new << std::endl;
    if (fraction_compatible_new < fraction_compatible) {
      double prob =
        std::exp((fraction_compatible_new - fraction_compatible) /
          factor_fraction_compatible);
      if (prob > rand(ctx.rng)) {
        if (ctx.get_verbose())
          std::cout << "Edge accepted despite of a lower fraction of "
                    << "compatible genotypes (accept. prob. = " << prob  << ")"
                    << std::endl;
        break;
      }
      if (i == p * p)
        std::cout << "All posible moves yielded a poset with lower fraction of"
                  << " compatible observations" << std::endl;
    } else {
      break;
    }
  }

  /* Use previous estimates as initialization for the EM */
  M_new.set_lambda(M.get_lambda());
  M_new.set_epsilon(M.get_epsilon());

  /* Compute likelihood of the proposed/new poset */
  llhood_new = MCEM_hcbn(
    M_new, obs, times, weights, L, sampling, version, control_EM,
    sampling_times_available, thrds, ctx);

  /* Accept the proposed poset, if it improves the likelihood */
  if (llhood_new > llhood_current) {
    num_accept += 1;
    M = M_new;
    fraction_compatible = fraction_compatible_new;
    if (ctx.get_verbose())
      std::cout << "Log-likelihood new poset:" << llhood_new << std::endl;
    return llhood_new;
  }
  
  /* Accept the proposed poset, if it doesn't improve the solution with certain
   * probability. However, when the temperature (T) is 0, only accept steps
   * that increase the likelihood
   */
  if (T == 0.0)
    return llhood_current;
  
  double acceptance_prob = std::exp(-(llhood_current - llhood_new)/T);
  
  if (acceptance_prob > rand(ctx.rng)) {
    num_accept += 1;
    M = M_new;
    fraction_compatible = fraction_compatible_new;
    if (ctx.get_verbose())
      std::cout << "Log-likelihood new poset:" << llhood_new << std::endl;
    return llhood_new;
  } else
    return llhood_current;
}

//' Simulated annealing
//' 
//' @noRd
//' @param lambdas rate parameters
double simulated_annealing(
    Model& poset, const MatrixXb& obs, const VectorXd& times,
    const RowVectorXd& weights, float T, const float adap_rate,
    const float acceptance_rate, const int step_size,
    const unsigned int max_iter_ASA, const float factor_fraction_compatible,
    const unsigned int L, const std::string& sampling,
    const unsigned int version, const ControlEM& control_EM,
    const std::string& output_dir, const bool sampling_times_available,
    const unsigned int thrds, Context& ctx) {

  const auto N = obs.rows();   // Number of observations / genotypes
  float acceptace_rate_current;
  float scaling_const = -std::log(2.0) / std::log(acceptance_rate);
  int num_accept = 0;

  /* 1. Compute likelihood of the initial model */
  double llhood = MCEM_hcbn(
    poset, obs, times, weights, L, sampling, version, control_EM,
    sampling_times_available, thrds, ctx);

  /* 2. Compute the fraction of compatible observations/genotypes with the
   *    initial poset
   */
  double fraction_compatible =
    (double) num_compatible_observations(obs, poset, N) / N;

  if (ctx.get_verbose()) {
    std::cout << "Initial Log-likelihood:" << llhood << std::endl;
    std::cout << "Fraction of compatible observations: " <<
      fraction_compatible << std::endl;
    std::cout << "Initial poset:\n";
    poset.print_cover_relations();
    std::cout << std::endl;
  }

  /* 3. (Adaptive) simulated annealing */
  // if (!boost::filesystem::is_directory(output_dir)) {
  //   std::cerr << "ERROR: Output directory '" << output_dir
  //             << "' does not exist!" << std::endl;
  //   return 0;
  // }
  std::ofstream outfile_temperature;
  outfile_temperature.open(output_dir + "asa.txt");
  outfile_temperature << "step\t llhood\t temperature\t acceptance rate"
                      << std::endl;

  std::ofstream outfile;
  outfile.open(output_dir + "params.txt");
  outfile << "step\t llhood\t epsilon\t lambdas" << std::endl;
  outfile << 0 << "\t" << llhood << "\t" << poset.get_epsilon() << "\t"
          << poset.get_lambda().transpose() << std::endl;

  for (unsigned int iter = 1; iter < max_iter_ASA; ++ iter) {
    if (ctx.get_verbose())
      std::cout << "Step " << iter << " - log-likelihood: " << llhood <<
        std::endl;

    /* 3.a Compute the likelihood for the new move */
    llhood = propose_edge(
      poset, obs, llhood, fraction_compatible, times, weights, T, num_accept,
      factor_fraction_compatible, L, sampling, version, control_EM,
      sampling_times_available, thrds, ctx);

    //DBG
    poset.print_cover_relations();
    /* Write to output file */
    outfile << iter << "\t" << llhood << "\t" << poset.get_epsilon() << "\t"
            << poset.get_lambda().transpose() << std::endl;

    /* 3.b Update temperature */
    if (!adaptive) {
      acceptace_rate_current = (float) num_accept / (iter + 1);
      T *= 1 - 1/poset.size();

      if (ctx.get_verbose())
        std::cout << "Temperature update: " << T << " (acceptance rate: "
                  << acceptace_rate_current << ")" << std::endl;
      /* Write to output file */
      outfile_temperature << iter << "\t" << llhood << "\t" << T << "\t"
                          << acceptace_rate_current << std::endl;
    } else if (!(iter % step_size)) {
      acceptace_rate_current = (float) num_accept / step_size;
      T *= std::exp((0.5 - std::pow(acceptace_rate_current, scaling_const)) *
        adap_rate);
      num_accept = 0;

      if (ctx.get_verbose())
        std::cout << "Temperature update: " << T << " (acceptance rate: "
                  << acceptace_rate_current << ")" << std::endl;
      /* Write to output file */
      outfile_temperature << iter << "\t" << llhood << "\t" << T << "\t"
                          << acceptace_rate_current << std::endl;
    }
  }

  /* 4. Compute likelihood of the final model */
  ControlEM control_last(control_EM.max_iter * 2,
                         control_EM.update_step_size * 2, tol, max_lambda);
  double llhood = MCEM_hcbn(
    poset, obs, times, weights, L, sampling, version, control_last,
    sampling_times_available, thrds, ctx);

  outfile_temperature.close();
  outfile.close();
  return llhood;
}

RcppExport SEXP _adaptive_simulated_annealing(
    SEXP posetSEXP, SEXP obsSEXP, SEXP lambdaSEXP, SEXP epsSEXP,
    SEXP timesSEXP, SEXP lambda_sSEXP, SEXP weightsSEXP, SEXP LSEXP,
    SEXP samplingSEXP, SEXP versionSEXP, SEXP max_iter_EMSEXP,
    SEXP update_step_sizeSEXP, SEXP tolSEXP, SEXP max_lambdaSEXP, SEXP T0SEXP,
    SEXP adap_rateSEXP, SEXP acceptance_rateSEXP, SEXP step_sizeSEXP,
    SEXP max_iter_ASASEXP, SEXP output_dirSEXP,
    SEXP sampling_times_availableSEXP, SEXP thrdsSEXP, SEXP verboseSEXP,
    SEXP seedSEXP) {
  
  try {
    /* Convert input to C++ types */
    // NOTE: Probably it is better to use Map, and return List with
    // poset and params
    MatrixXi poset = as<MapMati>(posetSEXP);
    // MatrixXi poset = as<MatrixXi>(posetSEXP);
    // MapMati poset(as<MapMati>(posetSEXP));
    const MatrixXb& obs = as<MatrixXb>(obsSEXP);
    const MapVecd lambda(as<MapVecd>(lambdaSEXP));
    const double eps = as<double>(epsSEXP);
    const MapVecd times(as<MapVecd>(timesSEXP));
    const float lambda_s = as<float>(lambda_sSEXP);
    const MapRowVecd weights(as<MapRowVecd>(weightsSEXP));
    const unsigned int L = as<unsigned int>(LSEXP);
    const std::string& sampling = as<std::string>(samplingSEXP);
    const unsigned int version = as<unsigned int>(versionSEXP);
    const unsigned int max_iter_EM = as<unsigned int>(max_iter_EMSEXP);
    const unsigned int update_step_size = as<unsigned int>(update_step_sizeSEXP);
    const double tol = as<double>(tolSEXP);
    const float max_lambda = as<float>(max_lambdaSEXP);
    const float T0 = as<float>(T0SEXP);
    const float adap_rate = as<float>(adap_rateSEXP);
    const float acceptance_rate = as<float>(acceptance_rateSEXP);
    const int step_size = as<int>(step_sizeSEXP);
    const unsigned int max_iter_ASA = as<unsigned int>(max_iter_ASASEXP);
    const std::string& output_dir = as<std::string>(output_dirSEXP);
    const bool sampling_times_available = as<bool>(sampling_times_availableSEXP);
    const int thrds = as<int>(thrdsSEXP);
    const bool verbose = as<bool>(verboseSEXP);
    const int seed = as<int>(seedSEXP);

    const auto p = poset.rows(); // Number of mutations / events
    const float factor_fraction_compatible = 0.05; //TODO: remove hard-coded value

    edge_container edge_list = adjacency_mat2list(poset);
    Model M(edge_list, p, lambda_s);
    M.set_lambda(lambda);
    M.set_epsilon(eps);
    M.has_cycles();
    if (M.cycle)
      throw not_acyclic_exception();
    M.topological_sort();

    ControlEM control_EM(max_iter_EM, update_step_size, tol, max_lambda);
    
    /* Call the underlying C++ function */
    Context ctx(seed, verbose);
    double llhood = adaptive_simulated_annealing(
      M, obs, times, weights, T0, adap_rate, acceptance_rate, step_size,
      max_iter_ASA, factor_fraction_compatible, L, sampling, version,
      control_EM, output_dir, sampling_times_available, thrds, ctx);
    
    /* Return the result as a SEXP */
    //TODO (return cover relations instead of adjacency matrix)
    return List::create(_["lambda"]=M.get_lambda(), _["eps"]=M.get_epsilon(),
                        _["poset"]=adjacency_list2mat(M), _["llhood"]=llhood);
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}

RcppExport SEXP _remove_test(SEXP posetSEXP) {

  try {
    /* Convert input to C++ types */
    MatrixXi poset = as<MapMati>(posetSEXP);

    const auto p = poset.rows(); // Number of mutations / events

    edge_container edge_list = adjacency_mat2list(poset);
    Model M(edge_list, p);
    M.has_cycles();
    if (M.cycle)
      throw not_acyclic_exception();
    M.topological_sort();

    /* Call the underlying C++ function */
    boost::graph_traits<Poset>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(M.poset); ei != ei_end; ++ei) {
      std::cout << "Testing edge: " << source(*ei, M.poset) << " --> "
                << target(*ei, M.poset) << std::endl;
      Node v1 = source(*ei, M.poset);
      Node v2 = target(*ei, M.poset);
      Model M_new(M);
      remove_edge(v1, v2, M_new.poset);
      M_new.has_cycles();
      if (M_new.cycle)
        throw not_acyclic_exception();
      M_new.topological_sort();
      if (!M_new.transitive_reduction_dag())
        std::cout << "Not transitively reduced" << std::endl;
    }

    /* Return the result as a SEXP */
    return wrap(0);
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}

RcppExport SEXP _transitive_reduction_dag(SEXP posetSEXP) {
  
  try {
    /* Convert input to C++ types */
    MatrixXi poset = as<MapMati>(posetSEXP);
    
    const auto p = poset.rows(); // Number of mutations / events
    
    edge_container edge_list = adjacency_mat2list(poset);
    Model M(edge_list, p);
    M.has_cycles();
    if (M.cycle)
      throw not_acyclic_exception();
    M.topological_sort();
    M.print_cover_relations();
    if (!M.transitive_reduction_dag()) {
      std::cout << "Not transitively reduced" << std::endl;
      M.print_cover_relations();
    }
    
    /* Return the result as a SEXP */
    return wrap(0);
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}
