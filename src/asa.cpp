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


void initialize_lambda(Model& model, const MatrixXb& obs, const float max_lambda) {

  unsigned int N = obs.rows();
  unsigned int p = model.size();
  VectorXd lambda(p);
  std::vector< std::unordered_set<Node> > all_predecessors(p);
  /* Loop through nodes in topological order */
  for (node_container::const_reverse_iterator v = model.topo_path.rbegin();
       v != model.topo_path.rend(); ++v) {
    /* Loop through (direct) predecessors/parents of node v */
    boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
    for (boost::tie(in_begin, in_end) = boost::in_edges(*v, model.poset);
         in_begin != in_end; ++in_begin) {
      Node u = source(*in_begin, model.poset);
      all_predecessors[*v].insert(u);
      all_predecessors[*v].insert(all_predecessors[u].begin(),
                                  all_predecessors[u].end());
    }

    /*float aux;
    int num_predecessors = all_predecessors[*v].size();
    if (num_predecessors > 0) {
      VectorXi partial_sums(N);
      partial_sums.setConstant(0);
      for (const auto& j: all_predecessors[*v])
        partial_sums = partial_sums + obs.col(j).cast<int>();
      aux = (partial_sums.array() == num_predecessors).count();
      partial_sums = partial_sums + obs.col(*v).cast<int>();
      aux = ((float) (partial_sums.array() == num_predecessors + 1).count()) / aux;
    } else { */
      /* Source node (i.e. without predecessors) */
      // TODO: Debug and benchmark
      /*aux = (float) obs.col(*v).sum() / N;
    } */
    /* Pseudocounts */
    unsigned int count = 1;
    float mutated = 1e-7;
    for (unsigned int i = 0; i < N; ++i) {
      bool parents_flag = true;
      for (const auto& j: all_predecessors[*v]) {
        if (!obs(i, j)) {
          parents_flag = false;
          break;
        }
      }
      if (parents_flag) {
        count += 1;
        if (obs(i, *v))
          mutated += 1.0;
      }
    }
    float aux = mutated / count;
    lambda[*v] = aux * model.get_lambda_s() / (1.0 - aux);
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
    const std::string& sampling, const unsigned int version,
    const ControlEM& control_EM, const bool sampling_times_available,
    const unsigned int thrds, Context& ctx) {
  
  const vertices_size_type p = M.size();  // Number of mutations / events
  const auto N = obs.rows();              // Number of observations / genotypes
  
  Model M_new(M);
  double llhood_new = 0.0;
  double fraction_compatible_new = 0.0;
  std::vector<int> idxs_pool(p);
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
  std::discrete_distribution<int> distribution({0.5, 0.1, 0.4});
  int move = distribution(ctx.rng);
  unsigned int i;

  switch(move) {
  case 0:
    /* Change edge
     * if the edge already existed, remove it. Otherwise, add it
     */
    for (i = 0; i < p * p; ++i) {

      /* Every time a new move is proposed, start with the initial poset */
      M_new = M;

      /* Draw one relation at random */
      unsigned int v1 = (int) cover_relations_pool[i] / p;
      unsigned int v2 = cover_relations_pool[i] % p;
      if (ctx.get_verbose())
        std::cout << "Testing edge: " << v1 << "->" << v2 << std::endl;

      if (v1 == v2)
        continue;

      std::pair<boost::graph_traits<Poset>::edge_descriptor, bool> add =
        boost::add_edge(v1, v2, M_new.poset);
      if (add.second) {
        /* Add edge - Check that the resulting graph doesn't contain cycles
         * and it is transitively reduced
         * NOTE: possible improvement - remove the edges that renders the
         * resulting graph not transitively reduced
         */
        M_new.has_cycles();
        if (M_new.cycle) {
          if (ctx.get_verbose())
            std::cout << "Cycle!" << std::endl;
          continue;
        } else {
          M_new.topological_sort();
          M_new.transitive_reduction_dag();
          if (M_new.reduction_flag) {
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
      if (heuristic_compatibility(M_new, obs, fraction_compatible,
                                  fraction_compatible_new,
                                  factor_fraction_compatible, ctx))
        break;
    }
    if (i == p * p)
      std::cout << "Warning: All posible moves yielded a poset with lower"
                << " fraction of compatible observations" << std::endl;
    break;
  case 1:
    /* Swap node labels: remove and add edge.
     * Pick random edge u -> v
     */
    for (i = 0; i <  boost::num_edges(M.poset) * 2; ++i) {

      /* Every time a new move is proposed, start with the initial poset */
      M_new = M;

      boost::graph_traits<Poset>::edge_descriptor e =
        boost::random_edge(M_new.poset, ctx.rng);
      Node u = source(e, M_new.poset);
      Node v = target(e, M_new.poset);

      if (ctx.get_verbose())
        std::cout << "Testing the edge swap: " << u << "->" << v << std::endl;
      /* Remove edge u -> v */
      remove_edge(u, v, M_new.poset);
      /* Add edge v -> u */
      std::pair<boost::graph_traits<Poset>::edge_descriptor, bool> add =
        boost::add_edge(v, u, M_new.poset);
      /* Check that the resulting graph doesn't contain cycles and it is
       * transitively reduced
       */
       M_new.has_cycles();
       if (M_new.cycle) {
         if (ctx.get_verbose())
           std::cout << "Cycle!" << std::endl;
         continue;
       } else {
         M_new.topological_sort();
         M_new.transitive_reduction_dag();
         if (ctx.get_verbose())
           std::cout << "Swapping edge: " <<  add.first << std::endl;
       }
       if (heuristic_compatibility(M_new, obs, fraction_compatible,
                                   fraction_compatible_new,
                                   factor_fraction_compatible, ctx))
         break;
    }
    if (i == boost::num_edges(M.poset) * 2)
      std::cout << "Warning: All tested swap moves yielded a poset with "
                << "lower fraction of compatible observations"
                << std::endl;
    break;
  case 2:
    /* Add or remove a cover relation. Add or remove an edge, while preserving
     * cover relations
     */
    for (i = 0; i < p * p; ++i) {
      M_new = M;

      /* Draw one relation at random */
      unsigned int v1 = (int) cover_relations_pool[i] / p;
      unsigned int v2 = cover_relations_pool[i] % p;
      if (ctx.get_verbose())
        std::cout << "Testing edge (preserving cover relations): " << v1 << "->"
                  << v2 << std::endl;

      if (v1 == v2)
        continue;

      M_new.set_children();
      const std::vector< std::unordered_set<Node> >& children = M_new.get_children();
      std::pair<boost::graph_traits<Poset>::edge_descriptor, bool> add =
        boost::add_edge(v1, v2, M_new.poset);
      if (add.second) {
        /* Add cover relation - Check that the resulting graph doesn't contain
        * cycles
        */
        M_new.has_cycles();
        if (M_new.cycle) {
          if (ctx.get_verbose())
            std::cout << "Cycle!" << std::endl;
          continue;
        } else {
          /* Check if cover relations before the edge addition are redundant */
          boost::graph_traits<Poset>::out_edge_iterator out_begin, out_end, next;
          boost::tie(out_begin, out_end) = out_edges(v1, M_new.poset);
          for (next = out_begin; out_begin != out_end; out_begin = next) {
            ++next;
            Node u = target(*out_begin, M_new.poset);
            if (children[v2].find(u) != children[v2].end()) {
              remove_edge(v1, u, M_new.poset);
              if (ctx.get_verbose())
                std::cout << "Removing redundant edge: " << v1 << "->" << u
                          << std::endl;
            }
          }
          if (ctx.get_verbose())
            M_new.print_cover_relations();

          M_new.topological_sort();
          M_new.transitive_reduction_dag();
          if (M_new.reduction_flag) {
            if (ctx.get_verbose())
              std::cout << "Adding new cover relation: " <<  add.first << std::endl;
          } else {
            continue;
          }
        }
      } else {
        /* Remove cover relation */
        remove_edge(v1, v2, M_new.poset);
        /* Update children of v1 */
        std::unordered_set<Node> v1_children = M_new.get_successors(v1);
        /* Check if edges of type v1 - > direct_children(v2) are required.
         * Possible improvement: check for edges pa(v1) to v2
         */
        boost::graph_traits<Poset>::out_edge_iterator out_begin, out_end;
        for (boost::tie(out_begin, out_end) = out_edges(v2, M_new.poset);
             out_begin != out_end; ++out_begin) {
          Node u = target(*out_begin, M_new.poset);
          if (v1_children.find(u) == v1_children.end()) {
            boost::add_edge(v1, u, M_new.poset);
            if (ctx.get_verbose())
              std::cout << "Adding lost cover relation: " << v1 << "->" << u
                        << std::endl;
          }
        }
        if (ctx.get_verbose())
          M_new.print_cover_relations();
        if (ctx.get_verbose())
          std::cout << "Removing cover relation: " <<  add.first << std::endl;
      }
      if (heuristic_compatibility(M_new, obs, fraction_compatible,
                                  fraction_compatible_new,
                                  factor_fraction_compatible, ctx))
        break;
    }
    break;
  }

  /* Initialization */
  initialize_lambda(M_new, obs, control_EM.max_lambda);
  M_new.set_epsilon((double) num_incompatible_events(obs, M_new) / (N * p));

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
    const RowVectorXd& weights, ControlSA& control_ASA, const unsigned int L,
    const std::string& sampling, const unsigned int version,
    const ControlEM& control_EM, const bool sampling_times_available,
    const unsigned int thrds, Context& ctx) {

  const auto N = obs.rows();   // Number of observations / genotypes
  float acceptace_rate_current;
  float scaling_const = -std::log(2.0) / std::log(control_ASA.get_acceptance_rate());
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
      version, control_EM, sampling_times_available, thrds, ctx);

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
  llhood = MCEM_hcbn(
    poset, obs, times, weights, L, sampling, version, control_last,
    sampling_times_available, thrds, ctx);

  outfile_temperature.close();
  outfile.close();
  return llhood;
}

RcppExport SEXP _adaptive_simulated_annealing(
    SEXP posetSEXP, SEXP obsSEXP, SEXP timesSEXP, SEXP lambda_sSEXP,
    SEXP weightsSEXP, SEXP LSEXP, SEXP samplingSEXP, SEXP versionSEXP,
    SEXP max_iter_EMSEXP, SEXP update_step_sizeSEXP, SEXP tolSEXP,
    SEXP max_lambdaSEXP, SEXP T0SEXP, SEXP adap_rateSEXP,
    SEXP acceptance_rateSEXP, SEXP step_sizeSEXP, SEXP max_iter_ASASEXP,
    SEXP adaptiveSEXP, SEXP outdirSEXP, SEXP sampling_times_availableSEXP,
    SEXP thrdsSEXP, SEXP verboseSEXP, SEXP seedSEXP) {
  
  try {
    /* Convert input to C++ types */
    // NOTE: Probably it is better to use Map, and return List with
    // poset and params
    MatrixXi poset = as<MapMati>(posetSEXP);
    // MatrixXi poset = as<MatrixXi>(posetSEXP);
    // MapMati poset(as<MapMati>(posetSEXP));
    const MatrixXb& obs = as<MatrixXb>(obsSEXP);
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

    ControlEM control_EM(max_iter_EM, update_step_size, tol, max_lambda);
    ControlSA control_ASA(outdir, acceptance_rate, T0, adap_rate, step_size,
                          max_iter_ASA, adaptive);

    /* Call the underlying C++ function */
    Context ctx(seed, verbose);
    double llhood = simulated_annealing(
      M, obs, times, weights, control_ASA, L, sampling, version, control_EM,
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
      remove_edge(v1, v2, M_new.poset);
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
