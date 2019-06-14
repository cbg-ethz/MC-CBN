/** mccbn: large-scale inference on conjunctive Bayesian networks
 *  Implementation of class Model
 *
 * @author Susana Posada Céspedes
 * @email susana.posada@bsse.ethz.ch
 */

#include <Rcpp.h>
#include <RcppEigen.h>
#include "mcem.h"
#include <vector>

unsigned int Model::size() const {
  return poset.rows();
}

void Model::fill_poset(const Eigen::Ref<const MatrixXi>& P) {
  poset = P;
}

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

VectorXd Model::get_lambda() const {
  return _lambda;
}

double Model::get_lambda(const unsigned int idx) const {
  return _lambda[idx];
}

float Model::get_lambda_s() const {
  return _lambda_s;
}

double Model::get_epsilon() const {
  return _epsilon;
}

double Model::get_llhood() const {
  return _llhood;
}

void Model::grandparents(unsigned int current, unsigned int parent,
                         std::vector<bool>& visit) {
  for (auto & k: pa_closure[parent]){
    if (!visit[k]) {
      pa_closure[current].push_back(k);
      N_pa_closure[current] += 1;
    }
    visit[k] = 1;
  }
}

void Model::grandchildren(unsigned int current, unsigned int child,
                          std::vector<bool>& visit) {
  for (auto & k: ch_closure[child]){
    if (!visit[k])
      ch_closure[current].push_back(k);
    visit[k] = 1;
  }
}

void Model::parents() {
  unsigned int p = N_pa.size();
  N_pa.setZero(p);
  
  // NOTE: the outer loop runs over the columns and the inner loop iterates
  // over the rows. This is done because Eigen stores matrices in column-major
  // order by default.
  for (unsigned int j = 0; j < p; ++j) {
    pa[j].clear();
    for (unsigned int i = 0; i < p; ++i) {
      if (poset(i, j) == 1) {
        N_pa(j) += 1;
        pa[j].push_back(i);
      }
    }
  }
}

void Model::parents(const Eigen::Ref<const VectorXi>& topological_path) {
  unsigned int p = N_pa.size();
  N_pa_closure = N_pa;
  pa_closure = pa;
  
  for (unsigned int i = 0; i < topological_path.size(); ++i) {
    unsigned int k = topological_path[i];
    std::vector<bool> visit(p);
    for (auto & j : pa[k])
      grandparents(k, j, visit);
  }
}

void Model::children(const Eigen::Ref<const VectorXi>& topological_path) {
  unsigned int p = N_pa.size();
  unsigned int k;
  
  // Add immediate children
  for (unsigned int j = 0; j < p; ++j)
    for (unsigned int i = 0; i < p; ++i)
      if (poset(i, j) == 1) 
        ch_closure[i].push_back(j);
      
      // Add grandchildren, great-grandchildren and so on ..
      for (int i = topological_path.size() - 1; i >= 0; --i) {
        k = topological_path[i];
        std::vector<bool> visit(p);
        std::vector<int> ch = ch_closure[k];
        for (auto & j : ch)
          grandchildren(k, j, visit);
      }
}
