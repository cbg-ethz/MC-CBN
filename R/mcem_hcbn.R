#' @title Monte Carlo Expectation Maximization
#' @export
#' 
#' @description parameter estimation for the hidden conjunctive Bayesian network
#' model (H-CBN) via importance sampling
#'
#' @param lambda a vector containing initial values for the rate parameters
#' @param poset a matrix containing the cover relations
#' @param obs a matrix containing observations or genotypes, where each row
#' corresponds to a genotype vector whose entries indicate whether an event has
#' been observed (1) or not (0)
#' @param lambda.s rate of the sampling process. Defaults to `1.0`
#' @param L number of samples to be drawn from the proposal in the E-step
#' @param eps an optional initial value of the error rate parameter
#' @param sampling sampling scheme. OPTIONS: 'naive' - generate occurrence
#' times according to current rate parameters, and, from them, generate the
#' genotypes; 'add-remove' - generate genotypes from observed genotypes using a
#' two-steps proposal. First, make genotypes compatible with the poset by
#' either adding or removing mutations. Second, perturb this version by adding
#' or removing a mutation while yielding a compatible observation; 'backward' -
#' generate a pool of compatible genotypes according to current rate parameters
#' and sample K observations proportional to the Hamming distance
#' @param times an optional vector containing times at which genotypes were
#' observed
#' @param weights an optional vector containing observation weights
#' @param version an optional argument indicating which version of the
#' `add-remove` sampling scheme to use. Defaults to `3`
#' @param perturb.prob probability of perturbing a genotype. Genotypes are
#' perturbed in order to learn the error rate, \code{epsilon}. A genotype is
#' perturbed, if after drawing a random number, this is smaller than
#' `perturb.prob`, otherwise `X_start = X`. This option is used if `sampling`
#' is set to `add-remove`. Defaults to `0.3`
#' @param max.iter the maximum number of EM iterations. Defaults to 100
#' iterations
#' @param max.lambda an optional upper bound on the value of the rate
#' parameters. Defaults to `1e6`
#' @param thrds number of threads for parallel execution
#' @param verbose an optional argument indicating whether to output logging
#' information
#' @param seed seed for reproducibility
MCEM.hcbn <- function(
  lambda, poset, obs, lambda.s=1.0, L, eps=NULL,
  sampling=c('naive', 'add-remove', 'backward'), times=NULL, weights=NULL,
  version=3L, perturb.prob=0.3, max.iter=100L, burn.in=0.8, max.lambda=1e6,
  thrds=1L, verbose=FALSE, seed=NULL) {
  
  sampling <- match.arg(sampling)
  N <- nrow(obs)
  if (!is.integer(poset))
    poset <- matrix(as.integer(poset), nrow=nrow(poset), ncol=ncol(poset))
  
  if (!is.integer(obs))
    obs <- matrix(as.integer(obs), nrow=N, ncol=ncol(obs))
  
  if (is.null(times)) {
    times <- numeric(N)
    sampling.times.available <- FALSE
  } else {
    sampling.times.available <- TRUE
  }
  if (is.null(weights))
    weights <- rep(1, N)
  
  if (is.null(seed))
    seed <- sample.int(3e4, 1)
  
  if (is.null(eps)) {
    set.seed(seed)
    eps <- runif(1, 0.01, 0.3)
  }
  .Call('_MCEM_hcbn', PACKAGE = 'mccbn', lambda, poset, obs, times,
        lambda.s, eps, weights, L, sampling, version, perturb.prob, max.iter,
        burn.in,  max.lambda, sampling.times.available, as.integer(thrds),
        verbose, as.integer(seed))
}

#' @title Importance sampling
#' @export
#' 
#' @description compute the sufficient statistics in expectation using
#' importance sampling
#' 
#' @param genotype a binary vector indicating whether an event has been
#' observed (`1`) or not (`0`)
#' @param L number of samples to be drawn from the proposal
#' @param poset a matrix containing the cover relations
#' @param lambda a vector of the rate parameters
#' @param eps error rate
#' @param time optional argument specifying the sampling time
#' @param sampling sampling scheme. OPTIONS: 'naive' - generate occurrence
#' times according to current rate parameters, and, from them, generate the
#' genotypes; 'add-remove' - generate genotypes from observed genotypes using a
#' two-steps proposal. First, make genotypes compatible with the poset by
#' either adding or removing mutations. Second, perturb this version by adding
#' or removing a mutation while yielding a compatible observation; 'backward' -
#' generate a pool of compatible genotypes according to current rate parameters
#' and sample K observations proportional to the Hamming distance
#' @param version an integer indicating which version of the `add-remove`
#' sampling scheme to use.
#' @param perturb.prob probability of perturbing a genotype
#' @param dist.pool Hamming distance between `genotype` and the genotype pool.
#' This option is used if `sampling` is set to `add-remove`
#' @param Tdiff.pool Expected time differences for the genotype pool. This
#' option is used if `sampling` is set to `add-remove`
#' @param lambda.s rate of the sampling process. Defaults to `1.0`
#' @param seed seed for reproducibility
importance.weight <- function(
  genotype, L, poset, lambda, eps, time=NULL,
  sampling=c('naive', 'add-remove', 'backward'), version, perturb.prob,
  dist.pool=integer(0), Tdiff.pool=matrix(0), lambda.s=1.0, seed=NULL) {
  
  sampling <- match.arg(sampling)
  if (!is.integer(genotype))
    genotype <- as.integer(genotype)
  
  if (!is.integer(poset))
    poset <- matrix(as.integer(poset), nrow=nrow(poset), ncol=ncol(poset))
  
  if (!is.integer(version))
    version <- as.integer(version)
  
  if (is.null(time)) {
    time <- 0
    sampling.times.available <- FALSE
  } else {
    sampling.times.available <- TRUE
  }
  if (!is.integer(dist.pool))
    dist.pool <- as.integer(dist.pool)
  
  if (is.null(Tdiff.pool))
    Tdiff.pool <- matrix(0)
  
  if (is.null(seed))
    seed <- sample.int(3e4, 1)
  
  .Call('_importance_weight', PACKAGE = 'mccbn', genotype, L, poset,
        lambda, eps, time, sampling, version, perturb.prob, dist.pool,
        Tdiff.pool, lambda.s, sampling.times.available, as.integer(seed))
}

#' @title Complete-data Log-Likelihood
#' @export
#' 
#' @description compute the complete-data log-likelihood or, equivalently, the
#' hidden log-likelihood
#' 
#' @param lambda a vector of the rate parameters
#' @param eps error rate
#' @param Tdiff a matrix of expected time differences
#' @param dist a vector of expected Hamming distances
#' @param W a vector containing observation weights
complete.loglikelihood <- function(lambda, eps, Tdiff, dist, W) {
  .Call('_complete_log_likelihood', PACKAGE = 'mccbn', lambda, eps, Tdiff, dist,
        W)
}