#ifndef MCEM_H
#define MCEM_H

#include <Rcpp.h>
#include <vector>

Rcpp::NumericMatrix generateMutationTimes(const Rcpp::NumericVector &lambda, 
                                          const Rcpp::IntegerMatrix &poset,
                                          Rcpp::IntegerMatrix &O,
                                          const std::vector<double> &times, 
                                          const Rcpp::IntegerVector &topo_path, 
                                          int nrOfSamples, float lambda_s,
                                          bool sampling_time_available);

#endif
