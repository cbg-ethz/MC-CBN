/** mccbn: large-scale inference on conjunctive Bayesian networks
 *  MCEM for the hidden conjuctive Bayesian network model
 *
 * @author Susana Posada Céspedes
 * @email susana.posada@bsse.ethz.ch
 */

#ifndef ADD_REMOVE_HPP
#define ADD_REMOVE_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include "mcem.hpp"
#include <vector>

VectorXd scale_path_to_mutation(const Model& model);

RowVectorXb draw_sample(const RowVectorXb& genotype, const Model& model,
                        const unsigned int move, const VectorXd& remove_weight,
                        const VectorXd& add_weight, double& q_choice,
                        const int idx_remove, const int idx_add,
                        bool compatible);

MatrixXd generate_mutation_times(const MatrixXb& obs, const Model& model,
                                 VectorXd& dens, VectorXd& sampling_time,
                                 Context::rng_type& rng,
                                 const bool sampling_times_available=false);

VectorXd cbn_density_log(const MatrixXd& time, const VectorXd& lambda);

#endif