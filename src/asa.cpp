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
#include <boost/graph/random.hpp>
// #include <boost/graph/filtered_graph.hpp>
// #include <boost/filesystem.hpp>

#include "mcem.hpp"
#include "asa.hpp"
#include "not_acyclic_exception.hpp"

using namespace Rcpp;

float ControlSA::get_adap_rate() const {
  return _adap_rate;
}

float ControlSA::get_acceptance_rate() const {
  return _acceptance_rate;
}

unsigned int ControlSA::get_step_size() const {
  return _step_size;
}

unsigned int ControlSA::get_max_iter() const {
  return _max_iter;
}

float ControlSA::get_compatible_fraction_factor() const {
  return _compatible_fraction_factor;
}

bool ControlSA::get_adaptive() const {
  return _adaptive;
}

const std::string& ControlSA::get_outdir() const{
  return _outdir;
}

void write_poset(const Poset& poset, const std::string& outdir) {
  std::ofstream outfile_poset;
  outfile_poset.open(outdir + "poset.txt", std::ofstream::trunc);
  boost::graph_traits<Poset>::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(poset); ei != ei_end; ++ei)
    outfile_poset << boost::get(
        boost::get(&Event::event_id, poset), source(*ei, poset))
    << "\t" << boost::get(
        boost::get(&Event::event_id, poset), target(*ei, poset))
    << std::endl;
  outfile_poset.close();
}

void initialize_lambda(Model& model, const MatrixXb& obs, const float max_lambda) {

  unsigned int N = obs.rows();
  unsigned int p = model.size();
  VectorXd lambda(p);
  auto id = boost::get(&Event::event_id, model.poset);
  std::vector< std::unordered_set<Node> > all_predecessors(p);
  /* Loop through nodes in topological order */
  for (node_container::const_reverse_iterator v = model.topo_path.rbegin();
       v != model.topo_path.rend(); ++v) {
    /* Loop through (direct) predecessors/parents of node v */
    boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
    for (boost::tie(in_begin, in_end) = boost::in_edges(*v, model.poset);
         in_begin != in_end; ++in_begin) {
      Node u = boost::get(id, source(*in_begin, model.poset));
      all_predecessors[model.poset[*v].event_id].insert(u);
      all_predecessors[model.poset[*v].event_id].insert(all_predecessors[u].begin(),
                                                        all_predecessors[u].end());
    }

    unsigned int count = 1;
    float mutated = 1e-7;
    int num_predecessors = all_predecessors[model.poset[*v].event_id].size();
    if (num_predecessors > 0) {
      VectorXi partial_sums(N);
      partial_sums.setConstant(0);
      for (const auto& j: all_predecessors[model.poset[*v].event_id])
        partial_sums = partial_sums + obs.col(j).cast<int>();
      count += (partial_sums.array() == num_predecessors).count();
      partial_sums = partial_sums + obs.col(model.poset[*v].event_id).cast<int>();
      mutated += (partial_sums.array() == num_predecessors + 1).count();
    } else {
      /* Source node (i.e. without predecessors) */
      count += N;
      mutated += obs.col(model.poset[*v].event_id).cast<int>().sum();
    }
    /* Pseudocounts */
    // unsigned int count = 1;
    // float mutated = 1e-7;
    // for (unsigned int i = 0; i < N; ++i) {
    //   bool parents_flag = true;
    //   for (const auto& j: all_predecessors[model.poset[*v].event_id]) {
    //     if (!obs(i, j)) {
    //       parents_flag = false;
    //       break;
    //     }
    //   }
    //   if (parents_flag) {
    //     count += 1;
    //     if (obs(i, model.poset[*v].event_id))
    //       mutated += 1.0;
    //   }
    // }
    float aux = mutated / count;
    lambda[model.poset[*v].event_id] = aux * model.get_lambda_s() / (1.0 - aux);
  }
  model.set_lambda(lambda, max_lambda);
}

//' Heuristic based on the fraction of compatible observations with the
//' proposed poset
bool heuristic_compatibility(Model& M, const MatrixXb& obs,
                             const double fraction_compatible,
                             double& fraction_compatible_new,
                             const float factor_fraction_compatible,
                             Context& ctx) {

  std::uniform_real_distribution<> rand(0.0, 1.0);

  const auto N = obs.rows();
  fraction_compatible_new =
    (double) num_compatible_observations(obs, M, N) / N;
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
        return true;
    }
  } else {
    return true;
  }
  return false;
}

//' Propose moves for the simulated annealing
//'
//' @noRd
double propose_edge(
    Model& M, const MatrixXb& obs, const double llhood_current,
    double& fraction_compatible, const VectorXd& times,
    const RowVectorXd& weights, const float T, int& num_accept,
    const float factor_fraction_compatible, const unsigned int L,
    const std::string& sampling, const ControlEM& control_EM,
    const bool sampling_times_available, MatrixXd& llhood_addRemove,
    MatrixXd& llhood_swap, MatrixXd& llhood_preserve, const unsigned int thrds,
    Context& ctx) {
  
  const vertices_size_type p = M.size();  // Number of mutations / events
  const auto N = obs.rows();              // Number of observations / genotypes
  
  Model M_new(M);
  double llhood_new = 0.0;
  double fraction_compatible_new = 0.0;
  unsigned int num_edges = boost::num_edges(M.poset);
  std::vector<int> idxs_pool(num_edges);
  std::vector<int> cover_relations_pool(p * p);

  /* Shuffle pool of indices. Use to draw edges/cover relations */
  std::iota(idxs_pool.begin(), idxs_pool.end(), 0);
  std::shuffle(idxs_pool.begin(), idxs_pool.end(), ctx.rng);
  std::iota(cover_relations_pool.begin(), cover_relations_pool.end(), 0);
  std::shuffle(cover_relations_pool.begin(), cover_relations_pool.end(), ctx.rng);

  std::uniform_real_distribution<> rand(0.0, 1.0);

  /* Propose an update.
   * Pick move: modify - add/delete - edge, swap edge or add/delete edge while
   * preserving cover relations
   */
  std::discrete_distribution<int> distribution({0.5, 0.0, 0.4, 0.1});
  int move = distribution(ctx.rng);
  unsigned int i, v1, v2;

  switch(move) {
  case 0:
    /* Change edge
     * if the edge already existed, remove it. Otherwise, add it
     */
    for (i = 0; i < p * p; ++i) {

      /* Every time a new move is proposed, start with the initial poset */
      M_new = M;

      /* Draw one relation at random */
      v1 = (int) cover_relations_pool[i] / p;
      v2 = cover_relations_pool[i] % p;
      if (v1 == v2)
        continue;

      if (ctx.get_verbose())
        std::cout << "Testing edge: " <<  M_new.poset[v1].event_id << "->"
                  << M_new.poset[v2].event_id << std::endl;

      bool add = M_new.add_edge(v1, v2);
      if (add) {
        if (M_new.cycle) {
          if (ctx.get_verbose())
            std::cout << "Cycle!" << std::endl;
          continue;
        } else {
          if (M_new.reduction_flag) {
            if (ctx.get_verbose())
              std::cout << "Adding new edge: (" <<  M_new.poset[v1].event_id
                        << "," << M_new.poset[v2].event_id << ")" << std::endl;
          } else {
            continue;
          }
        }
      } else {
        /* Remove edge*/
        M_new.remove_edge(v1, v2);
        if (ctx.get_verbose())
          std::cout << "Removing edge: (" <<  M_new.poset[v1].event_id << ","
                    << M_new.poset[v2].event_id << ")" << std::endl;
      }
      if (heuristic_compatibility(M_new, obs, fraction_compatible,
                                  fraction_compatible_new,
                                  factor_fraction_compatible, ctx))
        break;
    }
    if (i == p * p)
      std::cout << "Warning: All add/remove moves yielded a poset with lower "
                << "fraction of compatible observations" << std::endl;

    if (llhood_addRemove(v1, v2) == 0.0) {
      /* Initialization */
      initialize_lambda(M_new, obs, control_EM.max_lambda);
      M_new.set_epsilon((double) num_incompatible_events(obs, M_new) / (N * p));

      /* Compute likelihood of the proposed/new poset */
      llhood_new = MCEM_hcbn(
        M_new, obs, times, weights, L, sampling, control_EM,
        sampling_times_available, thrds, ctx);
      llhood_addRemove(v1, v2) = llhood_new;
    } else {
      llhood_new = llhood_addRemove(v1, v2);
    }
    break;
  case 1:
    /* Swap edge: remove and add edge.
     * Pick random edge v1 -> v2
     */
    for (i = 0; i < num_edges; ++i) {

      /* Every time a new move is proposed, start with the initial poset */
      M_new = M;

      boost::graph_traits<Poset>::edge_iterator ei;
      ei = boost::next(boost::edges(M_new.poset).first, idxs_pool[i]);
      Node v1 = source(*ei, M_new.poset);
      Node v2 = target(*ei, M_new.poset);

      if (ctx.get_verbose())
        std::cout << "Testing the edge swap: " << M_new.poset[v1].event_id
                  << "->" << M_new.poset[v2].event_id << std::endl;
      /* Remove edge v1 -> v2 */
      M_new.remove_edge(v1, v2);
      /* Add edge v2 -> v1 */
      M_new.add_edge(v2, v1);
      if (M_new.cycle) {
        if (ctx.get_verbose())
          std::cout << "Cycle!" << std::endl;
        continue;
      } else {
         if (ctx.get_verbose())
           std::cout << "Swapping edge: (" <<  M_new.poset[v1].event_id << ","
                     << M_new.poset[v2].event_id << ")" << std::endl;
      }
      if (heuristic_compatibility(M_new, obs, fraction_compatible,
                                 fraction_compatible_new,
                                 factor_fraction_compatible, ctx))
        break;
    }
    if (i == num_edges)
      std::cout << "Warning: All swap-edge moves yielded a poset with lower "
                << "fraction of compatible observations" << std::endl;

    /* Initialization */
    initialize_lambda(M_new, obs, control_EM.max_lambda);
    M_new.set_epsilon((double) num_incompatible_events(obs, M_new) / (N * p));

    /* Compute likelihood of the proposed/new poset */
    llhood_new = MCEM_hcbn(
      M_new, obs, times, weights, L, sampling, control_EM,
      sampling_times_available, thrds, ctx);
    break;
  case 2:
    /* Add or remove a cover relation. Add or remove an edge, while preserving
     * cover relations
     */
    for (i = 0; i < p * p; ++i) {
      /* Draw one relation at random */
      v1 = (int) cover_relations_pool[i] / p;
      v2 = cover_relations_pool[i] % p;
      if (v1 == v2)
        continue;

      if (ctx.get_verbose())
        std::cout << "Testing edge (preserving cover relations): "
                  << M.poset[v1].event_id << "->" << M.poset[v2].event_id
                  << std::endl;

      if (M.get_update_children())
        M.set_children();
      M_new = M;

      bool add = M_new.add_relation(v1, v2);
      if (add) {
        if (M_new.cycle) {
          if (ctx.get_verbose())
            std::cout << "Cycle!" << std::endl;
          continue;
        } else {
          if (M_new.reduction_flag) {
            /* NOTE: The resulting poset can be non transitively reduced, e.g.,
             * if the new edge doesn't correspond to a new cover relation (v2
             * was already reachable from v1)
             */
            if (ctx.get_verbose())
              std::cout << "Adding new cover relation: ("
                        << M_new.poset[v1].event_id << ","
                        << M_new.poset[v2].event_id << ")" << std::endl;
          } else {
            continue;
          }
        }
      } else {
        /* Remove cover relation */
        M_new.remove_relation(v1, v2);
        if (ctx.get_verbose())
          std::cout << "Removing cover relation: (" <<  M_new.poset[v1].event_id
                    << "," << M_new.poset[v2].event_id << ")" << std::endl;
      }
      if (heuristic_compatibility(M_new, obs, fraction_compatible,
                                  fraction_compatible_new,
                                  factor_fraction_compatible, ctx))
        break;
    }
    if (i == p * p)
      std::cout << "Warning: All moves yielded a poset with lower fraction "
                << "of compatible observations" << std::endl;

    if (llhood_preserve(v1, v2) == 0.0) {
      /* Initialization */
      initialize_lambda(M_new, obs, control_EM.max_lambda);
      M_new.set_epsilon((double) num_incompatible_events(obs, M_new) / (N * p));

      /* Compute likelihood of the proposed/new poset */
      llhood_new = MCEM_hcbn(
        M_new, obs, times, weights, L, sampling, control_EM,
        sampling_times_available, thrds, ctx);
      llhood_preserve(v1, v2) = llhood_new;
    } else {
      llhood_new = llhood_preserve(v1, v2);
    }
    break;
  case 3:
    /* Swap node labels */
    for (i = 0; i < p * p; ++i) {

      M_new = M;

      /* Draw two nodes at random */
      v1 = (int) cover_relations_pool[i] / p;
      v2 = cover_relations_pool[i] % p;
      if (v1 == v2)
        continue;

      if (ctx.get_verbose())
        std::cout << "Testing node swap: " << M_new.poset[v1].event_id << ", "
                  << M_new.poset[v2].event_id << std::endl;

      M_new.swap_node(v1, v2);
      if (ctx.get_verbose())
        M_new.print_cover_relations();

      if (heuristic_compatibility(M_new, obs, fraction_compatible,
                                  fraction_compatible_new,
                                  factor_fraction_compatible, ctx))
        break;
    }
    if (i == p * p)
      std::cout << "Warning: All swap-node moves yielded a poset with lower "
                << "fraction of compatible observations" << std::endl;

    if (llhood_swap(v1, v2) == 0.0) {
      /* Initialization */
      initialize_lambda(M_new, obs, control_EM.max_lambda);
      M_new.set_epsilon((double) num_incompatible_events(obs, M_new) / (N * p));

      /* Compute likelihood of the proposed/new poset */
      llhood_new = MCEM_hcbn(
        M_new, obs, times, weights, L, sampling, control_EM,
        sampling_times_available, thrds, ctx);
      llhood_swap(v1, v2) = llhood_swap(v2, v1) = llhood_new;
    } else {
      llhood_new = llhood_swap(v1, v2);
    }
    break;
  }

  /* Accept the proposed poset, if it improves the likelihood */
  if (llhood_new > llhood_current) {
    num_accept += 1;
    M = M_new;
    fraction_compatible = fraction_compatible_new;
    if (ctx.get_verbose())
      std::cout << "Log-likelihood new poset:" << llhood_new << std::endl;
    llhood_addRemove.setZero();
    llhood_swap.setZero();
    llhood_preserve.setZero();
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
    llhood_addRemove.setZero();
    llhood_swap.setZero();
    llhood_preserve.setZero();
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
    const RowVectorXd& weights, ControlSA& control_ASA, const unsigned int L,
    const std::string& sampling, const ControlEM& control_EM,
    const bool sampling_times_available, const unsigned int thrds,
    Context& ctx) {

  const auto N = obs.rows();   // Number of observations / genotypes
  float acceptace_rate_current;
  float scaling_const = -std::log(2.0) / std::log(control_ASA.get_acceptance_rate());
  int num_accept = 0;

  MatrixXd llhood_addRemove = MatrixXd::Zero(poset.size(), poset.size());
  MatrixXd llhood_swap = MatrixXd::Zero(poset.size(), poset.size());
  MatrixXd llhood_preserve = MatrixXd::Zero(poset.size(), poset.size());

  /* 1. Compute likelihood of the initial model */
  double llhood = MCEM_hcbn(
    poset, obs, times, weights, L, sampling, control_EM,
    sampling_times_available, thrds, ctx);
  Model poset_ML(poset);
  double llhood_ML = llhood;

  /* 2. Compute the fraction of compatible observations/genotypes with the
   *    initial poset
   */
  double fraction_compatible =
    (double) num_compatible_observations(obs, poset, N) / N;

  if (ctx.get_verbose()) {
    std::cout << "Initial observed Log-likelihood:" << llhood << std::endl;
    std::cout << "Fraction of compatible observations: " <<
      fraction_compatible << std::endl;
    std::cout << "Initial poset:\n";
    poset.print_cover_relations();
    std::cout << std::endl;
  }

  /* 3. (Adaptive) simulated annealing */
  // if (!boost::filesystem::is_directory(control_ASA.get_outdir())) {
  //   std::cerr << "ERROR: Output directory '" << control_ASA.get_outdir()
  //             << "' does not exist!" << std::endl;
  //   return 0;
  // }
  std::ofstream outfile_temperature;
  outfile_temperature.open(control_ASA.get_outdir() + "asa.txt");
  outfile_temperature << "step\t llhood\t temperature\t acceptance rate"
                      << std::endl;

  std::ofstream outfile;
  outfile.open(control_ASA.get_outdir() + "params.txt");
  outfile << "step\t llhood\t epsilon\t lambdas" << std::endl;
  outfile << 0 << "\t" << llhood << "\t" << poset.get_epsilon() << "\t"
          << poset.get_lambda().transpose() << std::endl;

  for (unsigned int iter = 1; iter < control_ASA.get_max_iter(); ++iter) {
    if (ctx.get_verbose())
      std::cout << "Step " << iter << " - log-likelihood: " << llhood
                << std::endl;

    /* 3.a Compute the likelihood for the new move */
    llhood = propose_edge(
      poset, obs, llhood, fraction_compatible, times, weights, control_ASA.T,
      num_accept, control_ASA.get_compatible_fraction_factor(), L, sampling,
      control_EM, sampling_times_available, llhood_addRemove, llhood_swap,
      llhood_preserve, thrds, ctx);
    if (llhood > llhood_ML) {
      llhood_ML = llhood;
      poset_ML = poset;
      write_poset(poset.poset, control_ASA.get_outdir());
    }

    if (ctx.get_verbose())
      poset.print_cover_relations();

    /* Write to output file */
    outfile << iter << "\t" << llhood << "\t" << poset.get_epsilon() << "\t"
            << poset.get_lambda().transpose() << std::endl;

    /* 3.b Update temperature */
    if (!control_ASA.get_adaptive()) {
      acceptace_rate_current = (float) num_accept / (iter + 1);
      control_ASA.T *= 1 - 1.0 / poset.size();

      if (ctx.get_verbose())
        std::cout << "Temperature update: " << control_ASA.T
                  << " (acceptance rate: " << acceptace_rate_current << ")"
                  << std::endl;
      /* Write to output file */
      outfile_temperature << iter << "\t" << llhood << "\t" << control_ASA.T
                          << "\t" << acceptace_rate_current << std::endl;
    } else if (!(iter % control_ASA.get_step_size())) {
      acceptace_rate_current = (float) num_accept / control_ASA.get_step_size();
      control_ASA.T *= std::exp((0.5 - std::pow(acceptace_rate_current, scaling_const)) *
        control_ASA.get_adap_rate());
      num_accept = 0;

      if (ctx.get_verbose())
        std::cout << "Temperature update: " << control_ASA.T
                  << " (acceptance rate: " << acceptace_rate_current << ")"
                  << std::endl;
      /* Write to output file */
      outfile_temperature << iter << "\t" << llhood << "\t" << control_ASA.T
                          << "\t" << acceptace_rate_current << std::endl;
    }
  }

  /* 4. Compute likelihood of the final model */
  ControlEM control_last(control_EM.max_iter * 2,
                         control_EM.update_step_size * 2, control_EM.tol,
                         control_EM.max_lambda);
  poset = poset_ML;
  llhood = MCEM_hcbn(
    poset, obs, times, weights, L, sampling, control_last,
    sampling_times_available, thrds, ctx);

  outfile_temperature.close();
  outfile.close();
  return llhood;
}

RcppExport SEXP _adaptive_simulated_annealing(
    SEXP posetSEXP, SEXP obsSEXP, SEXP timesSEXP, SEXP lambda_sSEXP,
    SEXP weightsSEXP, SEXP LSEXP, SEXP samplingSEXP, SEXP max_iter_EMSEXP,
    SEXP update_step_sizeSEXP, SEXP tolSEXP, SEXP max_lambdaSEXP, SEXP T0SEXP,
    SEXP adap_rateSEXP, SEXP acceptance_rateSEXP, SEXP step_sizeSEXP,
    SEXP max_iter_ASASEXP, SEXP neighborhood_distSEXP, SEXP adaptiveSEXP,
    SEXP outdirSEXP, SEXP sampling_times_availableSEXP, SEXP thrdsSEXP,
    SEXP verboseSEXP, SEXP seedSEXP) {
  
  try {
    /* Convert input to C++ types */
    MatrixXi poset = as<MapMati>(posetSEXP);
    const MatrixXb& obs = as<MatrixXb>(obsSEXP);
    const MapVecd times(as<MapVecd>(timesSEXP));
    const float lambda_s = as<float>(lambda_sSEXP);
    const MapRowVecd weights(as<MapRowVecd>(weightsSEXP));
    const unsigned int L = as<unsigned int>(LSEXP);
    const std::string& sampling = as<std::string>(samplingSEXP);
    const unsigned int max_iter_EM = as<unsigned int>(max_iter_EMSEXP);
    const unsigned int update_step_size = as<unsigned int>(update_step_sizeSEXP);
    const double tol = as<double>(tolSEXP);
    const float max_lambda = as<float>(max_lambdaSEXP);
    const float T0 = as<float>(T0SEXP);
    const float adap_rate = as<float>(adap_rateSEXP);
    const float acceptance_rate = as<float>(acceptance_rateSEXP);
    const int step_size = as<int>(step_sizeSEXP);
    const unsigned int max_iter_ASA = as<unsigned int>(max_iter_ASASEXP);
    const unsigned int neighborhood_dist = as<unsigned int>(neighborhood_distSEXP);
    const bool adaptive = as<bool>(adaptiveSEXP);
    const std::string& outdir = as<std::string>(outdirSEXP);
    const bool sampling_times_available = as<bool>(sampling_times_availableSEXP);
    const int thrds = as<int>(thrdsSEXP);
    const bool verbose = as<bool>(verboseSEXP);
    const int seed = as<int>(seedSEXP);

    const auto p = poset.rows(); // Number of mutations / events
    edge_container edge_list = adjacency_mat2list(poset);
    Model M(edge_list, p, lambda_s);
    M.has_cycles();
    if (M.cycle)
      throw not_acyclic_exception();
    M.topological_sort();

    /* Initialization */
    initialize_lambda(M, obs, max_lambda);
    M.set_epsilon((double) num_incompatible_events(obs, M) / (obs.rows() * p));

    ControlEM control_EM(max_iter_EM, update_step_size, tol, max_lambda,
                         neighborhood_dist);
    ControlSA control_ASA(outdir, acceptance_rate, T0, adap_rate, step_size,
                          max_iter_ASA, adaptive);

    /* Call the underlying C++ function */
    Context ctx(seed, verbose);
    double llhood = simulated_annealing(
      M, obs, times, weights, control_ASA, L, sampling, control_EM,
      sampling_times_available, thrds, ctx);
    
    /* Return the result as a SEXP */
    /* NOTE: (possible improvement) return cover relations instead of adjacency
     * matrix
     */
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
      boost::remove_edge(v1, v2, M_new.poset);
      M_new.has_cycles();
      if (M_new.cycle)
        throw not_acyclic_exception();
      M_new.topological_sort();
      M_new.transitive_reduction_dag();
      if (!M_new.reduction_flag)
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
    M.transitive_reduction_dag();
    if (!M.reduction_flag) {
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
