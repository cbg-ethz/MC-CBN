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
#' been observed (\code{1}) or not (\code{0})
#' @param lambda.s rate of the sampling process. Defaults to \code{1.0}
#' @param L number of samples to be drawn from the proposal in the E-step
#' @param eps an optional initial value of the error rate parameter
#' @param sampling sampling scheme. OPTIONS: \code{"forward"} - generate
#' occurrence times according to current rate parameters, and, from them,
#' generate the genotypes; \code{"add-remove"} - generate genotypes from
#' observed genotypes using a two-steps proposal. First, make genotypes
#' compatible with the poset by either adding or removing mutations. Second,
#' perturb this version by adding or removing a mutation while yielding a
#' compatible observation; \code{"rejection"} - generate a pool of compatible
#' genotypes according to current rate parameters and sample \code{K}
#' observations proportional to the Hamming distance
#' @param times an optional vector containing times at which genotypes were
#' observed
#' @param weights an optional vector containing observation weights
#' @param version an optional argument indicating which version of the
#' \code{"add-remove"} sampling scheme to use. Defaults to \code{3}
#' @param max.iter the maximum number of EM iterations. Defaults to \code{100}
#' iterations
#' @param update.step.size number of EM steps after which the number of
#' samples, \code{L}, is doubled. \code{L} is increased, if the difference in
#' the parameter estimates between such consecutive batches is greater than the
#' tolerance level, \code{tol}
#' @param tol convergence tolerance for the error rate and the rate parameters.
#' The EM runs until the difference between the average estimates in the last
#' two batches is smaller than tol, or until \code{max.iter} is reached.
#' @param max.lambda an optional upper bound on the value of the rate
#' parameters. Defaults to \code{1e6}
#' @param thrds number of threads for parallel execution
#' @param verbose an optional argument indicating whether to output logging
#' information
#' @param seed seed for reproducibility
MCEM.hcbn <- function(
  lambda, poset, obs, lambda.s=1.0, L, eps=NULL,
  sampling=c('forward', 'add-remove', 'rejection'), times=NULL, weights=NULL,
  version=3L, max.iter=100L, update.step.size=20L, tol=0.001, max.lambda=1e6,
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

  if (update.step.size > max.iter)
    update.step.size <- as.integer(max.iter / 5)

  if (is.null(seed))
    seed <- sample.int(3e4, 1)

  if (is.null(eps)) {
    set.seed(seed)
    eps <- runif(1, 0.01, 0.3)
  }
  .Call('_MCEM_hcbn', PACKAGE = 'mccbn', lambda, poset, obs, times,
        lambda.s, eps, weights, as.integer(L), sampling, version,
        as.integer(max.iter), as.integer(update.step.size), tol, max.lambda,
        sampling.times.available, as.integer(thrds), verbose, as.integer(seed))
}

#' @title Importance sampling
#' @export
#'
#' @description compute the sufficient statistics in expectation using
#' importance sampling
#'
#' @param genotype either a binary vector or a matrix containing observations
#' or genotypes. Each row of the matrix corresponds to a genotype vector whose
#' entries indicate whether an event has been observed (\code{1}) or not
#' (\code{0})
#' @param L number of samples to be drawn from the proposal
#' @param poset a matrix containing the cover relations
#' @param lambda a vector of the rate parameters
#' @param eps error rate
#' @param time optional argument specifying the sampling time
#' @param sampling sampling scheme. OPTIONS: \code{"forward"} - generate
#' occurrence times according to current rate parameters, and, from them,
#' generate the genotypes; \code{"add-remove"} - generate genotypes from
#' observed genotypes using a two-steps proposal. First, make genotypes
#' compatible with the poset by either adding or removing mutations. Second,
#' perturb this version by adding or removing a mutation while yielding a
#' compatible observation; \code{"rejection"} - generate a pool of compatible
#' genotypes according to current rate parameters and sample \code{K}
#' observations proportional to the Hamming distance
#' @param version an integer indicating which version of the
#' \code{"add-remove"} sampling scheme to use.
#' @param dist.pool Hamming distance between \code{genotype} and the genotype
#' pool. This option is used if \code{sampling} is set to \code{"rejection"}
#' and \code{genotype} corresponds to a vector containing a single genotype
#' @param Tdiff.pool Expected time differences for the genotype pool. This
#' option is used if \code{sampling} is set to \code{"rejection"} and
#' \code{genotype} corresponds to a vector containing a single genotype
#' @param lambda.s rate of the sampling process. Defaults to \code{1.0}
#' @param thrds number of threads for parallel execution. This option is used
#' if \code{genotype} corresponds to a matrix of genotypes
#' @param seed seed for reproducibility
importance.weight <- function(
  genotype, L, poset, lambda, eps, time=NULL,
  sampling=c('forward', 'add-remove', 'rejection'), version=NULL,
  dist.pool=integer(0), Tdiff.pool=matrix(0), lambda.s=1.0, thrds=1L,
  seed=NULL) {

  sampling <- match.arg(sampling)
  if (is.matrix(genotype))
    if (!is.integer(genotype))
      genotype <-
        matrix(as.integer(genotype), nrow=nrow(genotype), ncol=ncol(genotype))
  else if (!is.integer(genotype))
    genotype <- as.integer(genotype)

  if (!is.integer(poset))
    poset <- matrix(as.integer(poset), nrow=nrow(poset), ncol=ncol(poset))

  if (is.null(time)) {
    time <- 0
    if (is.matrix(genotype))
      time <- numeric(nrow(genotype))
    sampling.times.available <- FALSE
  } else {
    sampling.times.available <- TRUE
    if (is.matrix(genotype))
      if (length(sampling.times) != nrow(genotype))
        stop("A vector of length ",  nrow(genotype), " is expected")
  }
  if (is.null(version))
    version <- 0
  else if (!is.integer(version))
    version <- as.integer(version)

  if (!is.integer(dist.pool))
    dist.pool <- as.integer(dist.pool)

  if (is.null(Tdiff.pool))
    Tdiff.pool <- matrix(0)

  if (is.null(seed))
    seed <- sample.int(3e4, 1)

  if (is.matrix(genotype))
    .Call('_importance_weight', PACKAGE = 'mccbn', genotype, L, poset, lambda,
          eps, time, sampling, version, lambda.s, sampling.times.available,
          as.integer(thrds),  as.integer(seed))
  else
    .Call('_importance_weight_genotype', PACKAGE = 'mccbn', genotype, L, poset,
          lambda, eps, time, sampling, version, dist.pool, Tdiff.pool, lambda.s,
          sampling.times.available, as.integer(seed))
}

#' @title Complete-data Log-Likelihood
#' @export
#'
#' @description compute the complete-data log-likelihood or, equivalently, the
#' hidden log-likelihood
#'
#' @param lambda a vector of the rate parameters
#' @param eps error rate, \eqn{\epsilon}
#' @param Tdiff a matrix of expected time differences
#' @param dist a vector of expected Hamming distances
#' @param W sum of weighted observations
complete.loglikelihood <- function(lambda, eps, Tdiff, dist, W=NULL) {

  if (is.null(W))
    W <- length(Tdiff)

  .Call('_complete_log_likelihood', PACKAGE = 'mccbn', lambda, eps, Tdiff, dist,
        W)
}
