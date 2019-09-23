/** mccbn: large-scale inference on conjunctive Bayesian networks
 *  Implementation of class Model
 *
 * @author Susana Posada Céspedes
 * @email susana.posada@bsse.ethz.ch
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "mcem.hpp"
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/transitive_closure.hpp>
#include <vector>

void Model::set_lambda(const Eigen::Ref<const VectorXd>& lambda) {
  _lambda = lambda;
}

void Model::set_lambda(const Eigen::Ref<const VectorXd>& lambda, const float max_lambda) {
  _lambda = lambda;
  for (unsigned int j = 0; j < lambda.size(); ++j)
    if (lambda[j] > max_lambda)
      _lambda[j] = max_lambda;
}

void Model::set_epsilon(const double eps) {
  _epsilon = eps;
}

void Model::set_llhood(const double llhood) {
  _llhood = llhood;
}

void Model::has_cycles() {
  // Check for cycles using strongly connected components as proxy
  std::vector<vertices_size_type> component(_size);
  vertices_size_type num_scc = boost::strong_components(
    poset,
    boost::make_iterator_property_map(
      component.begin(), boost::get(boost::vertex_index, poset)));
  if (num_scc != _size)
    cycle = true;
}

void Model::topological_sort() {
  /* NOTE: topological_sort() writes its output in reverse topological order
   * (because it is more efficient to implement it that way). Suggestion: use
   * std::deque<int> (e.g. "std::deque<int> topo_path;") as its output data
   * structure, because it supports constant time insertion at the front, which
   * reverses the ordering. Moreover, topological_sort() requires one of two
   * things: (1) a color property map (e.g.
   * typedef adjacency_list< vecS, vecS, directedS, color_property<> > Graph;)
   * so that the algorithm can mark vertices to kepp track of its progress
   * through the graph, or (2) a mapping from vertices to integers so that the
   * algorithm can create its own color map with an array. When vertices are
   * integers, one could use 'identity_property_map' as the vertex index map.
   * std::deque<int> topo_path;
   * boost::topological_sort(poset, std::front_inserter(topo_path),
   *                         vertex_index_map(indentity_property_map()));
   */
  topo_path.clear();
  boost::topological_sort(poset, std::back_inserter(topo_path));
}

template <typename PropertyMap>
void Model::print_cover_relations(PropertyMap name) {
  boost::graph_traits<Poset>::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(poset); ei != ei_end; ++ei)
    std::cout << get(name, source(*ei, poset)) << " --> "
              << get(name, target(*ei, poset)) << std::endl;
}

//' @description Obtain (direct) predecessors/parents per node
std::vector<node_container> Model::get_direct_predecessors() const {

  std::vector<node_container> parents(_size);
  /* Loop through nodes in topological order */
  for (node_container::const_reverse_iterator v = topo_path.rbegin();
       v != topo_path.rend(); ++v) {
    /* Loop through (direct) predecessors/parents of node v */
    boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
    for (boost::tie(in_begin, in_end) = boost::in_edges(*v, poset);
         in_begin != in_end; ++in_begin) {
      Node u = source(*in_begin, poset);
      parents[*v].push_back(u);
    }
  }
  return parents;
}

//' @description Obtain all successors per node.
void Model::set_children() {

  /* Loop through vertices in reverse topological order */
  for (node_container::const_iterator it = topo_path.begin();
       it != topo_path.end(); ++it) {
    Node u = *it;
    _children[u].clear();
    /* Loop through (direct) successor/children of node u */
    boost::graph_traits<Poset>::out_edge_iterator out_begin, out_end;
    for (boost::tie(out_begin, out_end) = out_edges(u, poset);
         out_begin != out_end; ++out_begin) {
      Node v = target(*out_begin, poset);
      _children[u].insert(v);
      _children[u].insert(_children[v].begin(), _children[v].end());
    }
  }
}

//' @description Obtain successors of a given node.
std::unordered_set<Node> Model::get_successors(Node u) {
  std::unordered_set<Node> successors;

  /* Loop through (direct) successor/children of node u */
  boost::graph_traits<Poset>::out_edge_iterator out_begin, out_end;
  for (boost::tie(out_begin, out_end) = out_edges(u, poset);
       out_begin != out_end; ++out_begin) {
    Node v = target(*out_begin, poset);
    successors.insert(v);
    successors.insert(_children[v].begin(), _children[v].end());
  }
  return successors;
}

//' @description Obtain direct successors per node. Successsors are sorted by
//' the topological ordering
std::vector<node_container> Model::get_direct_successors(node_container& topo_order) {
  std::vector<node_container> successors(_size);
  boost::graph_traits<Poset>::vertex_iterator v_begin, v_end;
  for (boost::tie(v_begin, v_end) = boost::vertices(poset); v_begin != v_end;
       ++v_begin) {
    typedef boost::graph_traits<Poset>::adjacency_iterator adj_iter;
    std::pair<adj_iter, adj_iter> pr = boost::adjacent_vertices(*v_begin, poset);
    successors[*v_begin].assign(pr.first, pr.second);
    std::sort(successors[*v_begin].begin(), successors[*v_begin].end(),
              boost::bind(std::less<vertices_size_type>(),
                          boost::bind(
                            boost::detail::subscript(topo_order), _1),
                            boost::bind(
                              boost::detail::subscript(topo_order), _2)));
  }
  return successors;
}

//' @description Compute transitive reduction for a DAG
//' Code adapted from boost/graph/transitive_closure.hpp
void Model::transitive_reduction_dag() {
  reduction_flag = true;
  node_container topo_ordering(_size);
  vertices_size_type n = 0;
  /* Record the rank for each vertex within the topological ordering */
  for (node_container::reverse_iterator it = topo_path.rbegin();
       it != topo_path.rend(); ++it, ++n) {
    topo_ordering[*it] = n;
  }

  /* Direct successor per node and sorted according to the topological order */
  std::vector<node_container> successors = get_direct_successors(topo_ordering);

  /* Loop through vertices in reverse topological order and keep track of all
  * the successors of each vertex (direct and indirect successors). That is,
  * the set of all successors of a vertex v is the union of its direct and
  * indirect successors. The set of indirect successors of v is the union of
  * the successors of its direct successors.
  */
  std::vector< std::unordered_set<Node> > all_successors(_size);
  for (node_container::iterator it = topo_path.begin(); it != topo_path.end();
       ++it) {
    Node u = *it;
    node_container::const_iterator adj, adj_last;
    for (adj = successors[u].begin(), adj_last = successors[u].end();
         adj != adj_last; ++adj) {
      Node v = *adj;
      if (all_successors[u].find(v) != all_successors[u].end()) {
        /* remove edge (u, v) as v is already a children of u */
        remove_edge(u, v, poset);
        reduction_flag = false;
      } else {
        /* successors(u) = successors(u) U {v} */
        all_successors[u].insert(v);
        /* successors(u) = successors(u) U successors(v) */
        all_successors[u].insert(all_successors[v].begin(), all_successors[v].end());
      }
    }
  }
}

void Model::print_cover_relations() {
  print_cover_relations(boost::get(boost::vertex_index, poset));
}

void Model::clear() {
  topo_path.clear();
  cycle = false;
  reduction_flag = false;
}
