/** mccbn: large-scale inference on conjunctive Bayesian networks
 *
 * This file is part of the mccbn package
 * 
 * @author Susana Posada Céspedes
 * @email susana.posada@bsse.ethz.ch
 */

#include "mcem.hpp"
#include <boost/math/special_functions/binomial.hpp>

void draw_hidden_genotypes(MatrixXb& samples) {
  
  const unsigned int p = samples.rows();
  for (unsigned int j = 0; j < p; ++j)
    samples(j, j) = !samples(j, j);
  
}

unsigned int n_choose_k(unsigned int n, unsigned int k) {
  if (k > n)
    return 0;
  // if (k == 0 || k == n)
  //   return 1;
  // if (k * 2 > n)
  //   k = n - k;
  // 
  // unsigned int result = n;
  // for (unsigned int i = 2; i <= k; ++i) {
  //   result *= (n - i + 1);
  //   result /= i;
  // }
  // return result;
  return boost::math::binomial_coefficient<double>(n, k);
}

void neighbors(const unsigned int p, unsigned int k,
               Eigen::Ref<MatrixXb> output_mat) {

  if (p >= 1) {
    unsigned int nrows = n_choose_k(p - 1, k);
    unsigned int nrows_alt = n_choose_k(p - 1, k - 1);

    /* Mutation as observed */
    if (nrows > 0)
      neighbors(p - 1, k, output_mat.block(0, 1, nrows, p - 1));

    /* Change the observed mutation */
    if (nrows_alt > 0) {
      output_mat.block(nrows, 0, nrows_alt, 1).setConstant(!output_mat(nrows, 0));
      neighbors(p - 1, k - 1, output_mat.block(nrows, 1, nrows_alt, p - 1));
    }
  }
}

void runif_std(const unsigned int N, std::vector<double>& output,
               Context::rng_type& rng) {

  std::uniform_real_distribution<double> distribution(0, 1);
  for (unsigned int i = 0; i < N; ++i)
    output[i] = distribution(rng);
}

void coin_tossing(MatrixXb& output_mat, const double eps,
                  Context::rng_type& rng) {

  const unsigned int p = output_mat.cols();
  const unsigned int N = output_mat.rows();
  std::vector<double> result(p * N);
  runif_std(result.size(), result, rng);

  MatrixXd mask(p, N);
  mask = Map<MatrixXd>(result.data(), N, p);

  for (unsigned int i = 0; i < p; ++i)
    output_mat.col(i) = (mask.col(i).array() < eps).select(!output_mat(0, i), output_mat.col(i));
}

RcppExport SEXP _neighbors(SEXP pSEXP, SEXP kSEXP, SEXP genotypeSEXP) {
  try {
    const unsigned int p = Rcpp::as<unsigned int>(pSEXP);
    unsigned int k = Rcpp::as<unsigned int>(kSEXP);
    const RowVectorXb& genotype = Rcpp::as<RowVectorXb>(genotypeSEXP);
    
    std::vector<int> nrows(k);
    unsigned int nrows_sum = 0;
    for (unsigned int i = 0; i <= k; ++i) {
      nrows[i] = n_choose_k(p, i);
      nrows_sum += nrows[i];
    }
    MatrixXb samples = genotype.replicate(nrows_sum, 1);
    nrows_sum = 1;
    for (unsigned int i = 1; i <= k; ++i) {
      neighbors(p, i, samples.block(nrows_sum, 0, nrows[i], p));
      nrows_sum += nrows[i];
    }
    
    return Rcpp::wrap( samples );
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}

RcppExport SEXP _coin_tossing(SEXP epsSEXP, SEXP LSEXP, SEXP genotypeSEXP,
                              SEXP seedSEXP) {
  try {
    const double eps = Rcpp::as<double>(epsSEXP);
    const unsigned int L = Rcpp::as<unsigned int>(LSEXP);
    const RowVectorXb& genotype = Rcpp::as<RowVectorXb>(genotypeSEXP);
    const int seed = Rcpp::as<int>(seedSEXP);

    Context ctx(seed);
    MatrixXb samples = genotype.replicate(L, 1);
    coin_tossing(samples, eps, ctx.rng);

    return Rcpp::wrap( samples );
  } catch  (...) {
    handle_exceptions();
  }
  return R_NilValue;
}
