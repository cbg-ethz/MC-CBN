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
