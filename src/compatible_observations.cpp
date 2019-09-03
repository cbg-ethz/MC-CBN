/** mccbn: large-scale inference on conjunctive Bayesian networks
 *
 * This file is part of the mccbn package
 * 
 * @author Susana Posada Céspedes
 * @email susana.posada@bsse.ethz.ch
 */

#include "mcem.hpp"
#include "not_acyclic_exception.hpp"

bool is_compatible(const RowVectorXb& genotype, const Model& model) {

  /* Loop through nodes in reverse topological order */
  for (node_container::const_iterator v = model.topo_path.begin();
       v != model.topo_path.end(); ++v) {
    if (genotype[*v]) {
      /* Loop through (direct) predecessors/parents of node v */
      boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
      for (boost::tie(in_begin, in_end) = boost::in_edges(*v, model.poset);
           in_begin != in_end; ++in_begin)
        if (!genotype[source(*in_begin, model.poset)])
          return false;
    }
  }
  return true;
}

/*int num_incompatible(const RowVectorXb& genotype, const Model& poset) {

  int count = 0;
  boost::graph_traits<Poset>::vertex_iterator v_begin, v_end;
  for (boost::tie(v_begin, v_end) = boost::vertices(poset.poset);
       v_begin != v_end; ++v_begin) {
    Node v = *v_begin;
    boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
    for (boost::tie(in_begin, in_end) = boost::in_edges(*v_begin, poset.poset);
         in_begin != in_end; ++in_begin) {*/
      /* There is an edge u -> v. Genotype is not compatible with the poset if
       * event v is observed, while event u has not occurred yet
       */
      /*Node u = source(*in_begin, poset.poset);
      if (!genotype(u) && genotype(v))
        count += 1;
    }
  }
  return count;
}*/

int num_incompatible_events(const MatrixXb& genotype, const Model& poset) {

  int count = 0;
  boost::graph_traits<Poset>::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(poset.poset); ei != ei_end; ++ei) {
    /* There is an edge u -> v. Genotype is not compatible with the poset if
     * event v is observed, while event u has not occurred yet
     */
    Node u = source(*ei, poset.poset);
    Node v = target(*ei, poset.poset);
    VectorXi aux = genotype.col(u).cast<int>() - genotype.col(v).cast<int>();
    count += (aux.array() < 0).count();
  }
  return count;
}

int num_compatible_observations(const MatrixXb& obs, const Model& poset,
                                const unsigned int N) {
  
  int count = 0;
  for (unsigned int i = 0; i < N; ++i)
    count += is_compatible(obs.row(i), poset);

  return count;
}
  
  
std::vector<bool> compatible_observations(const MatrixXb& obs,
                                          const Model& poset,
                                          const unsigned int N) {
  
  std::vector<bool> idxs(N, false);
  for (unsigned int i = 0; i < N; ++i)
    if (is_compatible(obs.row(i), poset))
      idxs[i] = true;

  return idxs;
}

RcppExport SEXP _compatible_observations(SEXP obsSEXP, SEXP posetSEXP) {
  
  try {
    // Convert input to C++ types
    const MatrixXb& obs = Rcpp::as<MatrixXb>(obsSEXP);
    const MapMati poset(Rcpp::as<MapMati>(posetSEXP));
    
    const auto N = obs.rows();   // Number of observations / genotypes
    const auto p = poset.rows(); // Number of mutations / events
    edge_container edge_list = adjacency_mat2list(poset);
    Model M(edge_list, p);
    M.has_cycles();
    if (M.cycle)
      throw not_acyclic_exception();
    M.topological_sort();
    
    // Call the underlying C++ function
    std::vector<bool> ret = compatible_observations(obs, M, N);
    
    // Return the result as a SEXP
    return Rcpp::wrap( ret );
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}
