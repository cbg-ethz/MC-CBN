/** mccbn: large-scale inference on conjunctive Bayesian networks
 *
 * This file is part of the mccbn package
 * 
 * @author Susana Posada Céspedes
 * @email susana.posada@bsse.ethz.ch
 */

#include <boost/core/noncopyable.hpp>

#include "rng_utils.hpp"

#ifdef MKL_ENABLED

template <>
VectorXd rtexp<MklRng>(const unsigned int N, const double rate,
                                  const double cutoff, MklRng& rng) {
  std::vector<double> result(N), temp(N);
  
  VSLStreamStatePtr stream = rng.get_stream();
  
  if (VSL_STATUS_OK != vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE,
                                    stream, N, temp.data(), 0.0, 1.0))
    throw std::runtime_error("Something went wrong!");
  cblas_daxpy(N, std::expm1(-cutoff * rate), temp.data(), 1, result.data(), 1);
  
  std::swap(temp, result);
  vdLog1p(N, temp.data(), result.data());
  
  std::swap(temp, result);
  std::fill(result.begin(), result.end(), 0.0);
  cblas_daxpy(N, -1.0 / rate, temp.data(), 1, result.data(), 1);
  
  VectorXd output(result.size());
  output = Map<VectorXd>(result.data(), result.size());
  return output;
}

#else

template <>
VectorXd rtexp<StdRng>(const unsigned int N, const double rate,
                                  const double cutoff, StdRng& rng) {
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  VectorXd result(N);

  for (unsigned int i = 0; i < N; ++i) {
    double u = distribution(rng.get_rng());
    result[i] = -std::log1p(u * std::expm1(-cutoff * rate)) / rate;
  }

  return result;
}

#endif
