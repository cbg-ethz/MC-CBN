/** mccbn: large-scale inference on conjunctive Bayesian networks
 *  MCEM for the hidden conjuctive Bayesian network model
 *
 * @author Susana Posada Céspedes
 * @email susana.posada@bsse.ethz.ch
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <vector>
#include "mcem.hpp"
#include "not_acyclic_exception.hpp"

int rdiscrete(const VectorXd& weights, Context::rng_type& rng) {
  std::discrete_distribution<int> distribution(weights.data(),
                                               weights.data() + weights.size());
  return distribution(rng);
}

VectorXd scale_path_to_mutation(const Model& model) {
  
  VectorXd scale_cumulative(model.size());
  /* Loop through nodes in topological order */
  for (node_container::const_reverse_iterator v = model.topo_path.rbegin();
       v != model.topo_path.rend(); ++v) {
    scale_cumulative[*v] = 1 / model.get_lambda()[*v];
    double max_scale = -1;
    /* Loop through (direct) predecessors/parents of node v */
    boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
    for (boost::tie(in_begin, in_end) = boost::in_edges(*v, model.poset);
         in_begin != in_end; ++in_begin)
      max_scale =
        std::max(max_scale, scale_cumulative[source(*in_begin, model.poset)]);
      
    if (max_scale > -1)
      scale_cumulative[*v] += max_scale;
  }
  return scale_cumulative;
}

void add_all(RowVectorXb& genotype, const Model& model) {
  
  /* Loop through nodes in reverse topological order */
  for (node_container::const_iterator v = model.topo_path.begin();
       v != model.topo_path.end(); ++v) {
    if (genotype[*v]) {
      /* Loop through (direct) predecessors/parents of node v */
      boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
      for (boost::tie(in_begin, in_end) = boost::in_edges(*v, model.poset);
           in_begin != in_end; ++in_begin)
        if (!genotype[source(*in_begin, model.poset)])
          genotype[source(*in_begin, model.poset)] = true;
    }
  }
}

void add_parents(RowVectorXb& genotype, const Node v, const Model& model) {

  if (genotype[v]) {
    /* Loop through predecessors/parents of node v */
    boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
    for (boost::tie(in_begin, in_end) = boost::in_edges(v, model.poset);
         in_begin != in_end; ++in_begin) {
      Node u = source(*in_begin, model.poset);
      if (!genotype[u])
        genotype[u] = true;
      add_parents(genotype, u, model);
    }
  }
}

void remove_all(RowVectorXb& genotype, const Model& model) {

  /* Loop through nodes in topological order */
  for (node_container::const_reverse_iterator v = model.topo_path.rbegin();
       v != model.topo_path.rend(); ++v) {
    if (!genotype[*v]) {
      /* Loop through (direct) successor/children of node v */
      boost::graph_traits<Poset>::out_edge_iterator out_begin, out_end;
      for (boost::tie(out_begin, out_end) = out_edges(*v, model.poset);
           out_begin != out_end; ++out_begin)
        if (genotype[target(*out_begin, model.poset)])
          genotype[target(*out_begin, model.poset)] = false;
    }
  }
}

void remove_children(RowVectorXb& genotype, const Node v, const Model& model) {
  
  if (!genotype[v]) {
    /* Loop through (direct) successor/children of node v */
    boost::graph_traits<Poset>::out_edge_iterator out_begin, out_end;
    for (boost::tie(out_begin, out_end) = out_edges(v, model.poset);
         out_begin != out_end; ++out_begin) {
      Node u = target(*out_begin, model.poset);
      if (genotype[u])
        genotype[u] = false;
      remove_children(genotype, u, model);
    }
  }
}

RowVectorXb draw_sample(const RowVectorXb& genotype, const Model& model,
                        const bool move, const VectorXd& q_prob,
                        double& q_choice, Context::rng_type& rng) {

  int idx;
  VectorXd q_prob_inverse = q_prob.array().inverse();
  bool compatible = is_compatible(genotype, model);
  RowVectorXb sample = genotype;
  switch (move) {
    case 0 :
      /* Pick an event to be added */
      idx = rdiscrete(q_prob_inverse, rng);
      q_choice = q_prob_inverse[idx] / q_prob_inverse.sum();
      sample[idx] = true;
      if (compatible)
        add_parents(sample, idx, model);
      else
        add_all(sample, model);
      break;
    case 1 :
      /* Pick an event to be removed */
      idx = rdiscrete(q_prob, rng);
      q_choice = q_prob[idx] / q_prob.sum();
      sample[idx] = false;
      if (compatible)
        remove_children(sample, idx, model);
      else
        remove_all(sample, model);
      break;
  }
  return sample;
}


double rtexp(double time, double rate, std::mt19937& rng) {
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  double u = distribution(rng);

  double ret =
    time > 0 ? -std::log1p( u * std::expm1(-time * rate)) / rate : 0.0;
  return ret;
}

/* Probability density function of an exponential distribution in log scale */
double dexp_log(double time, double rate) {
  double ret = std::log(rate) - rate * time;
  return ret;
}

/* Cumulative distribution function of an exponential distribution in log scale */
double pexp_log(double time, double rate) {
  double ret = - std::expm1(-time * rate);
  ret = std::log(ret);
  return ret;
}

MatrixXd generate_mutation_times(
    const MatrixXb& obs, const Model& model, VectorXd& dens,
    VectorXd& sampling_time, Context::rng_type& rng, 
    const bool sampling_times_available=false) {

  unsigned int N = obs.rows();
  MatrixXd T_events;
  MatrixXd T_events_sum;
  T_events.setZero(N, obs.cols());
  T_events_sum.setZero(N, obs.cols());

  /* Generate sampling times sampling_time ~ Exp(lambda_{s}) */
  if (!sampling_times_available)
    sampling_time = rexp(N, model.get_lambda_s(), rng);
  
  /* Loop through nodes in topological order */
  for (node_container::const_reverse_iterator v = model.topo_path.rbegin();
       v != model.topo_path.rend(); ++v) {
    for (unsigned int i = 0; i < N; ++i) {
      double T_max = 0.0;
      /* Loop through parents of node v */
      boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
      for (boost::tie(in_begin, in_end) = boost::in_edges(*v, model.poset);
           in_begin != in_end; ++in_begin)
        if (T_events_sum(i, source(*in_begin, model.poset)) > T_max)
          T_max = T_events_sum(i, source(*in_begin, model.poset));

      if (obs(i, *v)) {
        /* Mutation is observed
         * Z ~ TExp(lambda, 0, sampling_time - time{max parents})
         */
        // TODO: Draw many random numbers at once.
        double time = rtexp(1, model.get_lambda()[*v], sampling_time[i] - T_max, rng)[0];
        dens[i] += dexp_log(time, model.get_lambda()[*v]) -
          pexp_log(sampling_time[i] - T_max, model.get_lambda()[*v]);
        T_events_sum(i, *v) =  T_max + time;
      } else {
        VectorXd time = rexp(1, model.get_lambda()[*v], rng) ;
        dens[i] += dexp_log(time[0], model.get_lambda()[*v]);
        T_events_sum(i, *v) = std::max(sampling_time[i], T_max) + time[0];
      }
      T_events(i, *v) = T_events_sum(i, *v) - T_max;   
    }
  }
  return T_events;
}


VectorXd cbn_density_log(const MatrixXd& time, const VectorXd& lambda) {
  
  unsigned int nrows = time.rows();
  VectorXd ret(nrows);
  ret.setConstant(lambda.array().log().sum());
  ret -= time * lambda;
  // for (unsigned int i = 0; i < nrows; ++i)
  //   ret[i] -= (time.row(i) * lambda);

  return ret;
}

RcppExport SEXP _scale_path_to_mutation(SEXP posetSEXP, SEXP lambdaSEXP) {
  try {
    // Convert input to C++ types
    const MapMati poset(Rcpp::as<MapMati>(posetSEXP));
    const MapVecd lambda(Rcpp::as<MapVecd>(lambdaSEXP));
    
    const auto p = poset.rows(); // Number of mutations / events
    edge_container edge_list = adjacency_mat2list(poset);
    Model M(edge_list, p);
    M.set_lambda(lambda);
    M.has_cycles();
    if (M.cycle)
      throw not_acyclic_exception();
    M.topological_sort();
    
    // Call the underlying C++ function
    VectorXd scale_cumulative = scale_path_to_mutation(M);
    
    // Return the result as a SEXP
    return Rcpp::wrap( scale_cumulative );
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}

RcppExport SEXP _draw_sample(
    SEXP genotypeSEXP, SEXP posetSEXP, SEXP moveSEXP, SEXP q_probSEXP,
    SEXP seedSEXP) {
  try {
    // Convert input to C++ types
    const RowVectorXb& genotype = Rcpp::as<RowVectorXb>(genotypeSEXP);
    const MapMati poset(Rcpp::as<MapMati>(posetSEXP));
    const bool move = Rcpp::as<bool>(moveSEXP);
    const MapVecd q_prob(Rcpp::as<MapVecd>(q_probSEXP));
    const int seed = Rcpp::as<int>(seedSEXP);

    const auto p = poset.rows(); // Number of mutations / events
    edge_container edge_list = adjacency_mat2list(poset);
    Model M(edge_list, p);
    M.has_cycles();
    if (M.cycle)
      throw not_acyclic_exception();
    M.topological_sort();
    
    Context ctx(seed);
    // Call the underlying C++ function
    double q_choice;
    RowVectorXb sample = draw_sample(genotype, M, move, q_prob, q_choice,
                                     ctx.rng);

    // Return the result as a SEXP
    return Rcpp::List::create(Rcpp::Named("sample")=sample,
                              Rcpp::Named("q_choice")=q_choice);
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}
  
RcppExport SEXP _generate_mutation_times(
    SEXP obsSEXP, SEXP posetSEXP, SEXP lambdaSEXP, SEXP timeSEXP,
    SEXP lambda_sSEXP, SEXP sampling_times_availableSEXP, SEXP seedSEXP) {
  try {
    // Convert input to C++ types
    const MatrixXb& obs = Rcpp::as<MatrixXb>(obsSEXP);
    const MapMati poset(Rcpp::as<MapMati>(posetSEXP));
    const MapVecd lambda(Rcpp::as<MapVecd>(lambdaSEXP));
    VectorXd time = Rcpp::as<MapVecd>(timeSEXP);
    const float lambda_s = Rcpp::as<float>(lambda_sSEXP);
    const bool sampling_times_available = Rcpp::as<bool>(sampling_times_availableSEXP);
    const int seed = Rcpp::as<int>(seedSEXP);
    
    const auto p = poset.rows(); // Number of mutations / events
    edge_container edge_list = adjacency_mat2list(poset);
    Model M(edge_list, p, lambda_s);
    M.set_lambda(lambda);
    M.has_cycles();
    if (M.cycle)
      throw not_acyclic_exception();
    M.topological_sort();
    
    Context ctx(seed);
    // Call the underlying C++ function
    VectorXd proposal_dens(obs.rows());
    MatrixXd Tdiff = generate_mutation_times(obs, M, proposal_dens, time,
                                             ctx.rng, sampling_times_available);

    // Return the result as a SEXP
    return Rcpp::List::create(Rcpp::Named("Tdiff")=Tdiff,
                              Rcpp::Named("proposal_density")=proposal_dens);
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}

RcppExport SEXP _cbn_density_log(SEXP timeSEXP, SEXP lambdaSEXP) {
  try {
    const MapMatd time(Rcpp::as<MapMatd>(timeSEXP));
    const MapVecd lambda(Rcpp::as<MapVecd>(lambdaSEXP));

    VectorXd ret = cbn_density_log(time, lambda);
    
    return Rcpp::wrap( ret );
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}

