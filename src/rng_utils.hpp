/** mccbn: large-scale inference on conjunctive Bayesian networks
 *
 * This file is part of the mccbn package
 * 
 * @author Susana Posada Céspedes
 * @email susana.posada@bsse.ethz.ch
 */

#ifndef RNG_UTILS_HPP
#define RNG_UTILS_HPP

#include <limits>
#include <RcppEigen.h>

using Eigen::Map;
using Eigen::VectorXd;

#define MKL_ENABLED

template <typename RNG_TYPE>
VectorXd rtexp(const unsigned int N, const double rate, const VectorXd& cutoff,
               RNG_TYPE& rng);

template <typename RNG_TYPE>
VectorXd rexp(const unsigned int N, const double rate, RNG_TYPE& rng) {
  double inf = std::numeric_limits<double>::infinity();
  VectorXd aux(N);
  aux.setConstant(inf);
  return rtexp(N, rate, aux, rng);
}

#ifdef MKL_ENABLED

#include "mkl_cblas.h"
#include "mkl_vsl.h"
#include "mkl_vml.h"

class MklRng: private boost::noncopyable {
public:
  MklRng(unsigned int seed) : _stream(NULL) {
    vslNewStream(&_stream, VSL_BRNG_MT19937, seed);
  }
  
  ~MklRng() {
    vslDeleteStream(&_stream);
  }
  
  VSLStreamStatePtr get_stream() const {
    return _stream;
  }
  
  typedef unsigned int result_type;
  
  result_type operator()() {
    unsigned int result;
    if (VSL_STATUS_OK != viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD,
                                            _stream, 1, &result))
      throw std::runtime_error("Something went wrong!");
    return result;
  }
  
  static constexpr result_type min() {
    return std::numeric_limits<result_type>::min();
  }

  static constexpr result_type max() {
    return std::numeric_limits<result_type>::max();
  }

private:
  VSLStreamStatePtr _stream;
};

typedef MklRng rng_type;

template <>
VectorXd rtexp<MklRng>(const unsigned int N, const double rate,
                       const VectorXd& cutoff, MklRng& rng);

#else

#include <random>
#include <algorithm>
#include <functional>
#include <array>

class StdRng: private boost::noncopyable {
public:
  StdRng(unsigned int seed) {
    std::ranlux24_base seeder_rng = std::ranlux24_base(seed);
    
    // adapted from https://stackoverflow.com/a/15509942
    std::array<int, std::mt19937::state_size> seed_data;
    std::generate_n(seed_data.data(), seed_data.size(), std::ref(seeder_rng));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    
    _rng.seed(seq);
  }
  
  typedef std::mt19937::result_type result_type;
  
  result_type operator()() {
    return _rng();
  }
  
  static constexpr result_type min() {
    return std::mt19937::min();
  }
  
  static constexpr result_type max() {
    return std::mt19937::max();
  }
  
  std::mt19937& get_rng() {
    return _rng;
  }
  
private:
  std::mt19937 _rng;
};

typedef StdRng rng_type;

template <>
VectorXd rtexp<StdRng>(const unsigned int N, const double rate,
                       const VectorXd& cutoff, StdRng& rng);

#endif

#endif
