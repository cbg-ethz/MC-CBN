/** 
 * Network learning using adaptive simulated annealing
 *
 * This file is part of the mccbn package
 * 
 * @author Susana Posada Céspedes
 * @email susana.posada@bsse.ethz.ch
 */

#ifndef ASA_HPP
#define ASA_HPP

#include <Rcpp.h>
#include "mcem.hpp"

/* Class containing customisable options for the simulated annealing 
 * algorithm
 */
class ControlSA {
public:
  float T;

  ControlSA(unsigned int p, const std::string& outdir, float T=50.0,
            float adap_rate=0.3, unsigned int step_size=20,
            unsigned int max_iter=1000, bool adaptive=true,
            float compatible_fraction_factor=0.05) :
    T(T), _outdir(outdir), _adap_rate(adap_rate), _step_size(step_size),
    _max_iter(max_iter), _adaptive(adaptive),
    _compatible_fraction_factor(compatible_fraction_factor) {
    _acceptance_rate = 1.0 / p;
  }
  
  ControlSA(const std::string& outdir, float acceptance_rate, float T=50.0,
            float adap_rate=0.3, unsigned int step_size=20,
            unsigned int max_iter=1000, bool adaptive=true,
            float compatible_fraction_factor=0.05) :
    T(T), _outdir(outdir), _acceptance_rate(acceptance_rate),
    _adap_rate(adap_rate), _step_size(step_size), _max_iter(max_iter),
    _adaptive(adaptive),
    _compatible_fraction_factor(compatible_fraction_factor) {}

  inline float get_adap_rate() const;

  inline float get_acceptance_rate() const;

  inline unsigned int get_step_size() const;

  inline unsigned int get_max_iter() const;

  inline float get_compatible_fraction_factor() const;
  
  inline bool get_adaptive() const;
  
  inline const std::string& get_outdir() const;

protected:
  std::string _outdir;
  float _acceptance_rate;
  float _adap_rate;
  unsigned int _step_size;
  unsigned int _max_iter;
  bool _adaptive;
  float _compatible_fraction_factor;
};

#endif
