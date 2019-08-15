#' @title Sample genotypes
#' @export
#'
#' @description Generate observations from a given poset and given rates
#'
#' @param N number of samples
#' @param poset a matrix containing the cover relations
#' @param lambda a vector of the rate parameters
#' @param sampling.times an optional vector of sampling times per observation
#' @param lambda.s rate of the sampling process. Defaults to \code{1.0}
#' @param seed seed for reproducibility
sample.genotypes <- function(N, poset, lambda, sampling.times=NULL,
                             lambda.s=1.0, seed=NULL) {
  
  p <- nrow(poset)
  if (!is.integer(poset))
    poset <- matrix(as.integer(poset), nrow=p, ncol=p)
  
  if (is.null(sampling.times)) {
    sampling.times <- numeric(N)
    sampling.times.available <- FALSE
  } else {
    sampling.times.available <- TRUE
    if (length(sampling.times) == 1)
      sampling.times <- rep(sampling.times, N)
    else if (length(sampling.times) != N)
      stop("A vector of length ", N, " is expected")
  }
  if (is.null(seed))
    seed <- sample.int(3e4, 1)
  
  .Call("_sample_genotypes", PACKAGE = 'mccbn', as.integer(N), poset, lambda,
        matrix(0, nrow=N, ncol=p), sampling.times, lambda.s,
        sampling.times.available, as.integer(seed))
}

#' @title Sample mutation times
#' @export
#'
#' @description Generate observation times from a given poset and given rates
#'
#' @inheritParams sample.genotypes
sample.times <- function(N, poset, lambda, seed=NULL) {

  p <- nrow(poset)
  if (!is.integer(poset))
    poset <- matrix(as.integer(poset), nrow=p, ncol=p)

  if (is.null(seed))
    seed <- sample.int(3e4, 1)

  .Call("_sample_times", PACKAGE = 'mccbn', as.integer(N), poset, lambda,
        matrix(0, nrow=N, ncol=p), as.integer(seed))
}

#' @title Generate genotypes
#' @export
#'
#' @description Generate observations from mutations times
#'
#' @inheritParams sample.genotypes
#'
#' @param mutation.times a matrix of dimension (\code{n} x \code{p}), where
#' \code{n} is the number of genotypes to be generated and \code{p} is the
#' number of events (or mutations)
#' @param sampling.time an optional argument specifying the sampling time
generate.genotypes <- function(mutation.times, poset, sampling.time=NULL,
                               lambda.s=1.0, seed=NULL) {

  N <- nrow(mutation.times)
  p <- nrow(poset)
  if (!is.integer(poset))
    poset <- matrix(as.integer(poset), nrow=p, ncol=p)

  if (is.null(sampling.time)) {
    sampling.time <- numeric(N)
    sampling.times.available <- FALSE
  } else {
    sampling.times.available <- TRUE
    if (length(sampling.time) == 1)
      sampling.time <- rep(sampling.time, N)
    else if (length(sampling.time) != N)
      stop("A vector of length ", N, " is expected")
  }
  if (is.null(seed))
    seed <- sample.int(3e4, 1)

  .Call("_generate_genotypes", PACKAGE = 'mccbn', mutation.times, poset,
        sampling.time, lambda.s, sampling.times.available, as.integer(seed))
}
