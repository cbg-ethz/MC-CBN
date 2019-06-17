/** mccbn: large-scale inference on conjunctive Bayesian networks
 *
 * @author Susana Posada Céspedes
 * @email susana.posada@bsse.ethz.ch
 */

#ifndef NOT_ACYCLIC_EXCEPTION_HPP
#define NOT_ACYCLIC_EXCEPTION_HPP

#include <exception>
#include <string>

/**
 * Exception thrown when a cycle is found.
 */
class not_acyclic_exception : public std::exception {
public:
  not_acyclic_exception(
    const char *error = "Poset does not correspond to an acyclic graph") {
    errorMessage = error;
  }
  
  const char *what() const noexcept {
    return errorMessage.c_str();
  }
  
private:
  std::string errorMessage;
};

#endif