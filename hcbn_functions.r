#' @title Hamming distance
#' @description compute the Hamming distance between (binary) genotypes, as the
#' \eqn{\sum_{i} |xi - yi|}
#'
#' @param x,y  binary vectors. Each entry in the vector indicates whether an
#' event has been observed (\code{1}) or not (\code{0})
hamming.dist <- function(x, y) {
  dist <- sum(abs(x - y))
  #dist <- sum(x != y)
  return(as.integer(dist))
}

#' @param p number of events/mutations
possible.genotypes <- function(p) {
  m <- 2^p
  genotypes <- matrix(0, nrow=m, ncol=p)
  for (i in 1:m) {
    genotypes[i, ] <- gen.genotype(p, m - i + 1)
  }
  return(genotypes)
}

#' @param p number of events/mutations
#' @param idx index identifying a given observation/genotype
gen.genotype <- function(p, idx) {
  g <- numeric(p)
  for(i in p:1) {
    g[i] <- idx %% 2
    idx <- idx %/% 2
  }
  return(g)
}

#' @description write poset to file in the format expected by H-CBN
#'
#' @param poset a matrix containing the cover relations
#' @param filename name of the output file
#' @param datadir output directory
write.poset <- function(poset, filename, datadir) {

  p <- nrow(poset)
  outfile <- file.path(datadir, paste(filename, ".poset", sep=""))
  #  write header
  write(p, outfile)

  relations <- which(poset == 1)
  i <- relations %% p
  j <- (relations %/% p) + 1
  write.table(cbind(i, j), outfile,row.names=FALSE, col.names=FALSE,
              append=TRUE)

  #  write footer
  write(0, outfile, append=TRUE)
}

#' @title Complete-data Log-Likelihood
#' @description compute the complete-data log-likelihood or, equivalently, the
#' hidden log-likelihood
#'
#' @param lambdas vector of the rate parameters
#' @param eps error rate
#' @param Tdiff matrix of expected time differences
#' @param dist vector of expected Hamming distances
complete.loglikelihood_ <- function(lambdas, eps, Tdiff, dist) {

  p <- length(lambdas)
  N <- length(dist)

  llhood <- N * sum(log(lambdas)) - sum(lambdas * t(Tdiff))
  if (eps == 0) {
    llhood <-
      ifelse(all(dist == 0), llhood,
             llhood + log(eps + .Machine$double.eps) * sum(dist) +
               log(1 - (eps + .Machine$double.eps)) * sum(p - dist))
  } else {
    llhood <- llhood + log(eps) * sum(dist) + log(1-eps) * sum(p - dist)
  }
  return(llhood)
}

#' @title Observed Log-Likelihood
#' @description compute the observed log-likelihood
#'
#' @param obs.events a matrix containing observations or genotypes: each row
#' correponds to a genotype vector whose entries indicate whether an event has
#' been observed (\code{1}) or not (\code{0})
#' @param poset a matrix containing the cover relations
#' @param lambdas vector of the rate parameters
#' @param lambda.s rate of the sampling process
#' @param eps error rate
#' @param sampling.times an optional vector of sampling times per observation
#' @param L number of samples to be drawn from the proposal. Defaults to
#' \code{1000}
#' @param genotype.pool an optional list containing sampled genotypes and
#' corresponding occurrence times
#' @param seed seed for reproducibility
obs.loglikelihood_ <- function(
  obs.events, poset, lambdas, lambda.s, eps, sampling.times=NULL, L=1000,
  sampling=c('forward', 'add-remove', 'rejection'), perturb.prob=0.8,
  version="3", genotype.pool=NULL, exact=FALSE, seed=NULL) {

  sampling <- match.arg(sampling)
  if (is.null(seed))
    seed <- sample.int(1e9, 1)

  set.seed(seed)
  seeds <- sample.int(3e4, N)

  # p: number of events / mutations
  # N: number of observations / genotypes
  if (exact) {
    ### **TODO** ### this is not correct
    # p <- length(lambdas)
    # N <- nrow(obs.events)
    # llhood <- 0
    # for (j in 1:p) {
    #   mutated <- sum(obs.events[, j] == 1)
    #   Pr_X_0 <- lambda.s/(lambda.s + lambdas[j])
    #   Pr_X_1 <- lambdas[j]/(lambda.s + lambdas[j])
    #   llhood <- llhood + mutated * log(eps * Pr_X_0 + (1-eps) * Pr_X_1)
    #   llhood <- llhood + (N - mutated) * log((1-eps) * Pr_X_0 + eps * Pr_X_1)
    # }
  } else {
    prob_Y <- foreach(i=1:N, .combine='c', .packages="mccbn") %dopar% {
      if (sampling == "rejection") {
        d_pool <-
          apply(genotype.pool$samples, 1, hamming.dist, y=obs.events[i, ])
        Tdiff_pool <- genotype.pool$Tdiff
      } else {
        d_pool <- NULL
        Tdiff_pool <- NULL
      }
      prob <-
        prob.importance.sampling(
          genotype=obs.events[i, ], L=L, poset=poset, lambdas=lambdas,
          lambda.s=lambda.s, eps=eps, sampling.time=sampling.times[i],
          sampling=sampling, perturb.prob=perturb.prob, version=version,
          dist.pool=d_pool, Tdiff.pool=Tdiff_pool, seed=seeds[i])
      return(prob)
    }

    llhood <- sum(log(prob_Y))
  }

  return(llhood)
}

#' @description empirical computation of the probability of a genotype
#' according to the model
#'
#' @param N number of samples
#' @param poset a matrix containing the cover relations
#' @param lambdas vector of the rate parameters
#' @param lambda.s rate of the sampling process
#' @param genotype a binary vector indicating whether an event has been
#' observed (\code{1}) or not (\code{0})
#' @param eps error rate, \eqn{\epsilon}
probY.empirical <- function(N, poset, lambdas, lambda.s, genotype, eps) {
  # TODO: uses sample.genotypes
  simGenotypes <-
    sample_genotypes(N, poset, sampling_param=lambda.s, lambdas=lambdas,
                     eps=eps)
  return(sum(apply(simGenotypes$obs_events, 1,
                   function(x, y) all(x == y), y=genotype)) / N)
}

#' @description compute the occurrence times expirically for a given genotype
#'
#' @inheritParams probY.empirical
tdiff.empirical <- function(N, poset, lambdas, lambda.s, genotype, eps) {

  simGenotypes <-
    sample_genotypes(N, poset, sampling_param=lambda.s, lambdas=lambdas,
                     eps=eps)
  idx <- which(apply(simGenotypes$obs_events, 1, function(x, y) all(x == y),
                    y=genotype))

  return(apply(simGenotypes$T_events[idx, ], 2, mean))
}

#' @description compute the average distance between genotype Y (subject to
#' noise) and N possible true genotypes, X
#'
#' @inheritParams probY.empirical
dist.empirical <- function(N, poset, lambdas, lambda.s, genotype, eps) {

  p <- ncol(poset)
  simGenotypes <-
    sample_genotypes(N, poset, sampling_param=lambda.s, lambdas=lambdas,
                     eps=eps)
  idx <-
    which(apply(simGenotypes$obs_events, 1, function(x, y) all(x == y),
                y=genotype))

  return(sum(
    apply(simGenotypes$hidden_genotypes[idx, ], 1, hamming.dist, y=genotype)) /
      length(idx))
}

#' @description identify the predecessors (parents) for each event
#'
#' @param poset a matrix containing the cover relations
#' @return list of parents per event
get.parents <- function(poset) {
  parents <- apply(poset, 2, function(x) which(x == 1))
  # Following step is need for the empty poset
  if (length(parents) == 0) {
    parents <- vector("list", length=nrow(poset))
  }
  return(parents)
}

#' @description identify the successors (children) for each event
#'
#' @param poset a matrix containing the cover relations
#' @return list of children per event
get.children <- function(poset) {
  children <- apply(poset, 1, function(x) which(x == 1))
  # Following step is need for the empty poset
  if (length(children) == 0) {
    children <- vector("list", length=nrow(poset))
  }
  return(children)
}

#' @description add events - look for parent nodes which have not occurred, but
#' any of their children had, and insert a 1 (i.e., replace `0` in the parent
#' node by `1`)
#'
#' @param genotype a binary vector indicating whether an event has been
#' observed (\code{1}) or not (\code{0})
#' @param children list of children per event. It should have been computed on
#' the transitively closed poset
add.relations <- function(genotype, children) {

  p <- length(genotype)
  return(sapply(1:p,
                function(i, obs, c)
                  ifelse(obs[i] == 0, ifelse(any(obs[c[[i]]] == 1), 1, obs[i]),
                         obs[i]), obs=genotype, c=children))
}

#' @description remove events - look for nodes which have occurred, but any of
#' their parents have not occurred, and remove it (i.e., replace \code{1} in
#' the child node by \code{0})
#'
#' @inheritParams add.relations
#' @param parents list of parents per node. It should have been computed on the
#' transitively closed poset
remove.relations <- function(genotype, parents) {

  p = length(genotype)
  return(sapply(1:p,
                function(i, obs, p)
                  ifelse(obs[i] == 1,
                         ifelse(any(obs[p[[i]]] == 0), 0, obs[i]), obs[i]),
                obs=genotype, p=parents))
}

#' @param genotype a binary vector indicating whether an event has been
#' observed (\code{1}) or not (\code{0})
#' @param parents list of parents per node. It should have been computed on the
#' transitively closed poset
#' @param children list of children per event. It should have been computed on
#' the transitively closed poset
get.possible.moves <- function(genotype, parents, children) {
  ### **TODO** ### vectorize get.possible.moves
  # Get the size of the Exit set - mutations that can happen next
  set.size.add <- 0
  idxs.add <- NA
  idx <- which(genotype == 0)
  if (length(idx) > 0) {
    idxs.add <- idx[sapply(idx, function(i, x, p) all(x[p[[i]]] == 1),
                           x=genotype, p=parents)]
    set.size.add <- length(idxs.add)
  }

  # Get the size of the set of mutations that happened last - can be removed
  # and the genotype remains compatible with the poset
  set.size.remove <- 0
  idxs.remove <- NA
  idx <- which(genotype == 1)
  if (length(idx) > 0) {
    idxs.remove <- idx[sapply(idx, function(i, x, c) all(x[c[[i]]] == 0),
                              x=genotype, c=children)]
    set.size.remove <- length(idxs.remove)
  }

  return(list("set.size"= c("add"=set.size.add, "remove"=set.size.remove),
              "idxs.add"=idxs.add, "idxs.remove"=idxs.remove))
}

#' @description perturb genotypes. Draw (uniformily) a random number \code{x}
#' between \code{0} and \code{1}. If \code{x} is less than \code{perturb.prob},
#' then perturb the genotype by adding or removing an event, but ensuring that
#' the resulting observation remains compatible with the poset
#'
#' @param events a matrix of genotypes, each row correponds to a vector
#' indicating whether an event has been observed (\code{1}) or not (\code{0})
#' @param parents list of parents per node. It should have been computed on the
#' transitively closed poset
#' @param children list of children per event. It should have been computed on
#' the transitively closed poset
#' @param perturb.prob probability of perturbing a genotype. Defaults to
#' \code{0.8}
#' @param compatible logical. If \code{TRUE} indicates that the observed
#' genotype, \code{Y}, is comptible with the current poset. Defaults to
#' \code{FALSE}
#' @param version indicates which version of the \code{"add-remove"} sampling
#' scheme to use
#' @return List with following named arguments: (1) \code{X} - perturbed
#' version of compatible genotypes, which remain compatible with the current
#' poset, (2) \code{option.set.size} size of the exit set or set of
#' observations which can be removed. \code{NA} is returned when the sample was
#' not perturbed
perturb <- function(events, parents, children, perturb.prob=0.8,
                    compatible=FALSE, version=c("1", "2", "3")) {

  version <- match.arg(version)

  L <- nrow(events)  # number of samples
  p <- ncol(events)  # number of events/mutations
  new.obs <- events
  option.set.size <- rep(NA, L)

  perturb <- ifelse(runif(L, 0, 1) < perturb.prob, 1, 0)
  idxs.perturb <- which(perturb == 1)

  # When observation Y was already compatible, all rows are equivalent
  if (compatible || var(as.vector(events)) == 0) {
    # Get indices of mutations that can be either added or removed
    idxs <- get.possible.moves(events[1, ], parents, children)
  }
  # Deal with special cases: wild type and resistant type
  if (all(events == 0)) {
    if (version == "1") {
      add <- runif(length(idxs.perturb), 0, 1) < 0.5
    } else if (version == "2" || version == "3") {
      add <- rep(TRUE, length(idxs.perturb))
    }
    mutation.idx <- sample(idxs$idxs.add, sum(add), replace=TRUE)
    option.set.size[idxs.perturb[add]] <- idxs$set.size["add"]
    if (sum(add) == 1) {
      new.obs[idxs.perturb[add], mutation.idx] = 1
    } else {
      ### **TODO** ### vectorize it
      for (i in 1:sum(add)) {
        new.obs[idxs.perturb[add][i], mutation.idx[i]] = 1
      }
    }
  } else if (all(events == 1)) {
    if (version == "1") {
      add <- runif(length(idxs.perturb), 0, 1) < 0.5
    } else if (version == "2" || version == "3") {
      add <- rep(FALSE, length(idxs.perturb))
    }
    mutation.idx <- sample(idxs$idxs.remove, sum(!add), replace=TRUE)
    option.set.size[idxs.perturb[!add]] <- idxs$set.size["remove"]
    if (sum(!add) == 1) {
      new.obs[idxs.perturb[!add], mutation.idx] <- 0
    } else {
      ### **TODO** ### vectorize it
      for (i in 1:sum(!add)) {
        new.obs[idxs.perturb[!add][i], mutation.idx[i]] <- 0
      }
    }
  } else {
    for (i in idxs.perturb) {

      # Get indices of mutations that can be either added or removed
      if (!compatible) {
        idxs <- get.possible.moves(events[i, ], parents, children)
      }

      if (version == "1") {
        # Add or remove one observation, both moves have equal probability.
        # Draw a random number between 0 and 1. If this number is > 0.5, then
        # add an observation
        add <- runif(1, 0, 1) < 0.5
        if (add) {
          mutation.idx <-
            ifelse(idxs$set.size["add"] > 1, sample(idxs$idxs.add, 1), idxs$idxs.add)
          option.set.size[i] <- idxs$set.size["add"]
          new.obs[i, mutation.idx] <- 1
        } else {
          mutation.idx <-
            ifelse(idxs$set.size["remove"] > 1, sample(idxs$idxs.remove, 1),
                   idxs$idxs.remove)
          option.set.size[i] <- idxs$set.size["remove"]
          new.obs[i, mutation.idx] <- 0
        }
      } else if (version == "2" || version == "3") {
        option.set.size[i] <- sum(idxs$set.size)

        # Choose one index randomly
        idx <- sample.int(option.set.size[i], 1)
        mutation.idx <- na.omit(c(idxs$idxs.add, idxs$idxs.remove))[idx]

        # Check whether move corresponds to an adding move
        add <- ifelse(idx <= idxs$set.size["add"], TRUE, FALSE)

        new.obs[i, mutation.idx] <- ifelse(add, 1, 0)
      }
    }
  }
  return(list("X"=new.obs, "option.set.size"=option.set.size))
}

#' @description draw L samples from the proposal. Two-steps proposal: (1) make
#' compatible and (2) perturb genotype
#'
#' @param genotype a binary vector indicating whether an event has been
#' observed (\code{1}) or not (\code{0})
#' @param L number of samples to be drawn from the proposal
#' @param parents list of parents per node. It should have been computed on the
#' transitively closed poset
#' @param children list of children per event. It should have been computed on
#' the transitively closed poset
#' @param compatible variable indicating whether the genotype is compatible
#' (\code{1}) with current poset or not (\code{0})
#' @param eps error rate
#' @param perturb.prob probability of perturbing a genotype. Genotypes are
#' perturbed in order to learn the error rate, \eqn{\epsilon}. Defaults to
#' \code{0.8}
#' @param version indicates which version of the \code{"add-remove"} sampling
#' scheme to use
draw.samples <- function(genotype, L, parents, children, compatible, eps,
                         perturb.prob=0.8, version=c("1", "2", "3")) {

  version <- match.arg(version)

  p <- length(genotype)  # number of events/mutations

  if (compatible) {
    # if genotype is compatible with poset, we can inmediately go to step 2:
    # perturb
    compatible.obs <- matrix(genotype, nrow=L, ncol=p, byrow=TRUE)
    make.compatible.prob <- NA
  } else {
    # if genotype is NOT compatible with poset, we need to generate L compatible
    # versions by adding or removing observations
    compatible.by.adding <- add.relations(genotype, children)
    compatible.by.removing <- remove.relations(genotype, parents)

    if (version == "1" | version == "2") {
      # Both moves, add and remove, have equal probability
      # Move add (remove) wouldn't be applicable if genotype is the resistant
      # type (wild type). However, such genotypes are always compatible with
      # the poset
      add <- ifelse(runif(L, 0, 1) < 0.5, 1, 0)
      make.compatible.prob <- rep(0.5, L)
    } else if (version == "3") {
      # Choose to add or remove observations according to the number of
      # modifications
      dist.add <- hamming.dist(genotype, compatible.by.adding)
      dist.remove <- hamming.dist(genotype, compatible.by.removing)
      add.prob <- eps^dist.add / (eps^dist.add + eps^dist.remove)
      add <- ifelse(runif(L, 0, 1) < add.prob, 1, 0)
      make.compatible.prob <- ifelse(add == 1, add.prob, 1 - add.prob)
    }
    compatible.obs <- matrix(0, L, p)

    if (all(add == 1)) {
      compatible.obs <- matrix(compatible.by.adding, nrow=L, ncol=p, byrow=TRUE)
    } else if (all(add == 0)) {
      compatible.obs <-
        matrix(compatible.by.removing, nrow=L, ncol=p, byrow=TRUE)
    } else {
      idxs.add <- which(add == 1)
      compatible.obs[idxs.add, ] = matrix(compatible.by.adding, nrow=sum(add),
                                          ncol=p, byrow=TRUE)
      compatible.obs[-idxs.add, ] = matrix(compatible.by.removing, nrow=L-sum(add),
                                           ncol=p, byrow=TRUE)
    }
  }
  samples <-
    perturb(compatible.obs, parents, children, perturb.prob, compatible,
            version=version)
  return(list("samples"=samples, "make.compatible.prob"=make.compatible.prob))
}

#' @description random generation for the truncated exponential distribution
#'
#' @param x upper bound on the observation
#' @param rate rate of the exponential distribution (i.e. mean 1/`rate`)
rtexp <- function(x, rate) {
  rand_num <- runif(1, 0, 1)
  Z <- ifelse(x > 0, -log(1 - rand_num * (1 - exp(-rate * x))) / rate, 0)
  return(Z)
}

#' @param genotype a binary vector indicating whether an event has been
#' observed (\code{1}) or not (\code{0})
#' @param poset a matrix containing the cover relations
#' @param lambdas a vector of the rate parameters
#' @param lambda.s rate of the sampling process. Defaults to \code{1.0}
#' @param sampling.time optional argument specifying the sampling time
sample.mutation.times <- function(genotype, poset, lambdas, lambda.s=1.0,
                                  sampling.time=NULL) {

  p <- length(genotype)  # number of events/mutations
  parents <- get.parents(poset)
  if (is.null(sampling.time))
    sampling.time = rexp(1, rate=lambda.s)

  Tdiff <- numeric(p)
  Tsum <- numeric(p)
  dens <- 0
  topo.path <- my.topological.sort(poset)

  for (j in topo.path) {
    max.time.parents <- 0
    if (length(parents[[j]]) > 0) {
      for (u in 1:length(parents[[j]])) {
        if(Tsum[parents[[j]][u]] > max.time.parents)
          max.time.parents <- Tsum[parents[[j]][u]]
      }
    }

    if (genotype[j] == 1) {
      # If mutation is observed,
      # Z ~ TExp(lambda, 0, sampling.time - max.time.parents)
      # When a lot of mutations have occurred max.time.parents approaches the
      # 'sampling.time'. Due to numerical precision problems,
      # 'sampling.time' < 'max.time.parents'. Therefore, we have to make sure
      # that when it happens rtexp returns 0
      # dexp: P(Z)
      # pexp: P(Z <= sampling.time - max.time.parents)
      Z <- rtexp(sampling.time - max.time.parents, rate=lambdas[j])
      Tsum[j] <-  max.time.parents + Z
      dens <- dens +  dexp(Z, rate=lambdas[j], log=TRUE) -
        pexp(sampling.time - max.time.parents, rate=lambdas[j], log.p=TRUE)
    } else {
      # If mutation is not observed, Z ~ Exp(lambda)
      Z <- rexp(1, rate=lambdas[j])
      Tsum[j] <- max(sampling.time, max.time.parents) + Z
      dens <-  dens + dexp(Z, rate=lambdas[j], log=TRUE)
    }
    Tdiff[j] <- Tsum[j] - max.time.parents
  }
  return(c("density"=dens, "Tdiff"=Tdiff))
}

#' @param Tdiff matrix of expected time differences
#' @param rate vector of exponential rates
log.cbn.density <- function(Tdiff, rate) {
  dens <- sum(log(rate)) - sum(rate * Tdiff)
  return(dens)
}

#' @param genotype a binary vector indicating whether an event has been
#' observed (\code{1}) or not (\code{0})
#' @param L number of samples to be drawn from the proposal
#' @param poset a matrix containing the cover relations
#' @param lambdas a vector of exponential rates
#' @param lambda.s rate of the sampling process
#' @param eps error rate
#' @param sampling.time an optional argument specifying the sampling time
#' @param sampling type of sampling scheme. OPTIONS: \code{"forward"},
#' \code{"add-remove"} or \code{"rejection"}
#' @param perturb.prob probability of perturbing a genotype. Genotypes are
#' perturbed in order to learn the error rate, \eqn{\epsilon}. Defaults to
#' \code{0.8}
#' @param version an optional argument indicating which version of the
#' \code{"add-remove"} sampling scheme to use. Defaults to \code{3}
#' @param dist.pool Hamming distance between \code{genotype} and the genotype
#' pool. This option is used if \code{sampling} is set to \code{"add-remove"}
#' @param seed seed for reproducibility
#' @return returns a list containing the importance weights and the sufficient
#' statistics. If \code{sampling} is set to \code{"rejection"}, it returns the
#' importance weights and the indices of choosed observations from the genotype
#' pool
importance.weight_ <- function(
  genotype, L, poset, lambdas, lambda.s, eps, sampling.time=NULL,
  sampling=c('forward', 'add-remove', 'rejection'), perturb.prob=0.8,
  version="3", dist.pool=NULL, seed=NULL) {

  sampling <- match.arg(sampling)
  if (!is.null(seed))
    set.seed(seed)

  p <- ncol(poset) # number of events/mutations

  if (sampling == 'forward') {

    # Generate L samples from poset with parameters 'lambdas' and 'lambda.s'
    # In particular, epsilon is zero (default value) - because the idea is to
    # generate samples of X (underlying true)
    ### **TODO** ### what to do when sampling times available. At the moment,
    # not considered to generate samples - incorporate sampling times in
    # function sample_genotypes
    if (!is.null(sampling.time)) {
      warning("Forward proposal doesn't account for sampling times")
    }
    samples <- sample_genotypes(L, poset, sampling_param=lambda.s, lambdas=lambdas)
    d <- apply(samples$hidden_genotypes, 1, hamming.dist, y=genotype)
    prob.Y.X <- eps^d * (1-eps)^(p-d)

    return(list("w"=prob.Y.X, "time.differences"=samples$T_events, "dist"=d))

  } else if (sampling == 'add-remove') {

    poset.trans.closed <- trans_closure(poset)
    parents <- get.parents(poset.trans.closed)
    children <- get.children(poset.trans.closed)
    compatible <- is_compatible(genotype, poset)
    # Generate L samples according to proposal
    return.list <-
      draw.samples(genotype, L, parents, children, compatible, eps,
                   perturb.prob, version)
    make.compatible.prob <- return.list$make.compatible.prob
    samples <- return.list$samples

    # Generate mutation times Tj from sample i
    Tdiff <- apply(samples$X, 1, sample.mutation.times, poset=poset,
                   lambdas=lambdas, lambda.s=lambda.s,
                   sampling.time=sampling.time)

    log.proposal.X <- Tdiff["density", ]
    Tdiff <- t(tail(Tdiff, n=-1))

    # Hamming distance bewteen samples and genotype
    dist <- apply(samples$X, 1, hamming.dist, y=genotype)

    # Computing log(Pr(Y|X))
    if (eps == 0) {
      # NOTE: If all observations are compatible with poset and pertub_prob == 0
      #       (which can happen because noisy observations can be compatible),
      #       then eps can be 0 and dist 0
      log.Pr.Y.X <-
        ifelse(dist == 0, 0,
               dist * log(eps + .Machine$double.eps) +
                 (p - dist) * log(1 - (eps + .Machine$double.eps)))
    } else {
      log.Pr.Y.X <- dist*log(eps) + (p - dist) * log(1 - eps)
    }

    # Computing log(Pr(X))
    log.Pr.X <- apply(Tdiff, 1, log.cbn.density, rate=lambdas)

    # Computing density of the proposal - for correction
    # proposal : weight accounting for making genotype compatible + weight
    #            accounting for choosing a possible mutation for perturbing
    #            the sample + weight accounting for sampling times from a
    #            truncated exponential
    num.options <- samples$option.set.size
    if (version == "1") {
      log.proposal.Y.X <-
        ifelse(is.na(num.options) | num.options == 0, log(1 - perturb.prob),
               log(perturb.prob) + log(0.5) + log(1/num.options))
    } else if (version == "2" || version == "3") {
      stopifnot(all(na.omit(num.options) != 0))
      log.proposal.Y.X <- ifelse(is.na(num.options), log(1 - perturb.prob),
                                 log(perturb.prob) + log(1 / num.options))
    }

    if (!compatible)
      log.proposal.Y.X <- log.proposal.Y.X + log(make.compatible.prob)

    log.proposal <- log.proposal.Y.X + log.proposal.X

    importance.weight <- exp(log.Pr.X + log.Pr.Y.X - log.proposal)
    return(list("w"=importance.weight, "time.differences"=Tdiff, "dist"=dist))

  } else if (sampling == 'rejection') {
    ### **TODO** ### what to do when sampling times available. At the moment,
    # not considered to generate samples - incorporate sampling times in
    # function sample_genotypes
    if (is.null(dist.pool)) {
      stop("Vector of distances expected for rejection sampling")
    }
    # if (is.null(genotype_pool)) {
    #   stop("Pool of genotypes are expected for rejection sampling")
    # }
    K <- length(dist.pool)
    # compatible = is_compatible(genotype, poset)
    aux <- eps^(dist.pool) * (1-eps) ^ (p - dist.pool)
    # if (compatible) {
    #   genotype_pool <- rbind(genotype_pool, genotype)
    #   q.prob <- c(aux, 1 - sum(aux))
    #   idxs <- sample(seq(1, K + 1), L, replace=TRUE, prob=q.prob)
    # } else {
    #   q.prob <- aux/sum(aux)
    #   idxs <- sample(seq(1, K), L, replace=TRUE, prob=q.prob)
    # }
    q.prob <- aux / sum(aux)
    idxs <- sample(seq(1, K), L, replace=TRUE, prob=q.prob)
    # Generate mutation times Tj from sample i
    # Tdiff <-
    #   apply(genotype_pool[idxs, ], 1, sample.mutation.times, poset=poset,
    #         lambdas=lambdas, lambda.s=lambda.s, sampling.time=sampling.time)
    #
    # log.proposal.X = Tdiff["density", ]
    # Tdiff = t(tail(Tdiff, n=-1))

    # Hamming distance bewteen samples and genotype
    dist <- dist.pool[idxs]

    # Computing log(Pr(Y|X))
    if (eps == 0) {
      # NOTE: If all obs ervations are compatible with poset and pertub_prob == 0
      #       (which can happen because noisy observations can be compatible),
      #       then eps can be 0 and dist 0
      log.Pr.Y.X <-
        ifelse(dist == 0, 0,
               dist * log(eps + .Machine$double.eps) +
                 (p - dist) * log(1 - (eps + .Machine$double.eps)))
    } else {
      log.Pr.Y.X <- dist * log(eps) + (p - dist) * log(1 - eps)
    }

    # Computing log(Pr(X))
    # log.Pr.X <- apply(Tdiff, 1, log.cbn.density, rate=lambdas)

    # Computing density of the proposal - for correction
    log.proposal <- log(q.prob[idxs]) + log(K) #+ log.proposal.X

    # importance.weight <- exp(log.Pr.X + log.Pr.Y.X - log.proposal)
    importance.weight <- exp(log.Pr.Y.X - log.proposal)
    # return(list("w"=importance.weight, "time.differences"=Tdiff, "dist"=dist))
    return(list("w"=importance.weight, "idxs"=idxs))
  }
}

#' @description compute Pr(Y) using Monte Carlo sampling, where `Y` corresponds
#' to the observed genotype
#'
#' @inheritParams importance.weight_
#' @param Tdiff.pool Expected time differences for the genotype pool. This
#' option is used if \code{sampling} is set to \code{"add-remove"}
prob.importance.sampling <- function(
  genotype, L, poset, lambdas, lambda.s, eps, sampling.time=NULL,
  sampling=c('forward', 'add-remove', 'rejection'), perturb.prob=0.8,
  version="3", dist.pool=NULL, Tdiff.pool=NULL, seed=NULL) {

  sampling <- match.arg(sampling)
  if (is.null(seed))
    seed <- sample.int(3e4, 1)

  if (sampling == "add-remove") {
    set.seed(seed, kind="L'Ecuyer-CMRG")
    probs <-
      importance.weight_(genotype, L, poset, lambdas, lambda.s, eps,
                         sampling.time=sampling.time, sampling=sampling,
                         perturb.prob=perturb.prob, version=version,
                         dist.pool=dist.pool)
  } else {
    if (is.null(Tdiff.pool))
      Tdiff.pool <- matrix(0)
    probs <-
      importance.weight(genotype, L, poset, lambdas, eps, time=sampling.time,
                        sampling=sampling, version=version, dist.pool=dist.pool,
                        Tdiff.pool=Tdiff.pool, lambda.s=lambda.s,
                        seed=as.integer(seed))
  }
  probs <- probs$w

  return(sum(probs) / L)
}

#' @description compute expected time differences for a given genotype using
#' Monte Carlo sampling
#'
#' @inheritParams importance.weight_
tdiff.importance.sampling <- function(
  genotype, L, poset, lambdas, lambda.s, eps, sampling.time=NULL,
  sampling=c('forward', 'add-remove', 'rejection'), perturb.prob=0.8,
  version="3", dist.pool=NULL, seed=NULL) {

  sampling <- match.arg(sampling)
  if (!is.null(seed))
    set.seed(seed)

  importance.weights <-
    importance.weight_(genotype, L, poset, lambdas, lambda.s, eps,
                       sampling.time=sampling.time, sampling=sampling,
                       perturb.prob=perturb.prob, version=version,
                       dist.pool=dist.pool)
  return(colSums(importance.weights$w * importance.weights$time.differences) /
           sum(importance.weights$w))

}

#' @description compute expected Hamming distance between a given observation
#' (genotype) and the underlying/true genotype using Monte Carlo sampling
#'
#' @inheritParams importance.weight_
dist.importance.sampling <- function(
  genotype, L, poset, lambdas, lambda.s, eps, sampling.time=NULL,
  sampling=c('forward', 'add-remove', 'rejection'), perturb.prob=0.8,
  version="3", dist.pool=NULL, seed=NULL) {

  sampling <- match.arg(sampling)
  if (!is.null(seed))
    set.seed(seed)

  importance.weights <-
    importance.weight_(genotype, L, poset, lambdas, lambda.s, eps,
                       sampling.time=sampling.time, sampling=sampling,
                       perturb.prob=perturb.prob, version=version,
                       dist.pool=dist.pool)
  return(sum(importance.weights$dist * importance.weights$w) /
           sum(importance.weights$w))
}

#' @param events a matrix of observations or a vector containing a single
#' genotype, if \code{one.genotype} is set to \code{TRUE}
#' @param L number of samples to be drawn from the proposal
#' @param poset a matrix containing the cover relations
#' @param lambdas a vector of exponential rates
#' @param lambda.s rate of the sampling process
#' @param eps error rate
#' @param rep an optional argument indicating the number of repetitions. This
#' option is used if \code{one.genotype} is set to \code{TRUE}
#' @param one.genotype a boolean variable indicating whether only one genotype
#' is provided
#' @param sampling.times an optional vector of sampling times per observation
#' @param sampling type of sampling scheme. OPTIONS: \code{"forward"},
#' \code{"add-remove"} or \code{"rejection"}
#' @param perturb.prob probability of perturbing a genotype. This option is
#' used if \code{sampling} is set to \code{"add-remove"}. Defaults to
#' \code{0.8}
#' @param version an optional agrument indicating which version of the
#' \code{"add-remove"} sampling scheme to use. Defaults to \code{3}
#' @param genotype.pool an optional list containing sampled genotypes and their
#' corresponding occurrence times
#' @param outdir an optional argument indicating the path to the output
#' directory
#' @param outname an optional argument indicating the suffix to include in the
#' output file name
#' @param binwidth an optional argument indicating the width of the histrogram
#' bins. This option is used if \code{one.genotype} is set to \code{TRUE}
#' @param seed seed for reproducibility
prob.empirical_vs_sampling <- function(
  events, L, poset, lambdas, lambda.s, eps, rep=NULL, one.genotype=FALSE,
  sampling.times=NULL, sampling=c('forward', 'add-remove', 'rejection'),
  perturb.prob=0.8, version="3", genotype.pool=NULL, outdir=NULL, outname="",
  binwidth=0.01, seed=NULL) {

  sampling <- match.arg(sampling)
  if (is.null(seed))
    seed <- sample.int(1e9, 1)

  set.seed(seed)
  seeds <- sample.int(3e4, N)

  if (one.genotype) {
    if (length(sampling.times) > 1) {
      warning("Only one sampling time was expected. First entry of vector ",
               "\'sampling.times\' is used.")
      sampling.times <- sampling.times[1]
    }
    genotype <- events
    sampling.time <- sampling.times
    N <- rep
  } else {
    N <- nrow(events)
  }

  if (sampling == "rejection") {
    if (one.genotype)
      d_pool <- apply(genotype.pool$samples, 1, hamming.dist, y=genotype)
    Tdiff_pool <- genotype.pool$Tdiff
  } else {
    d_pool <- NULL
    Tdiff_pool <- NULL
  }

  probs <- foreach(i=1:N, .combine='rbind', .packages="mccbn") %dopar% {
    if (!one.genotype) {
      genotype <- events[i, ]
      sampling.time <- sampling.times[i]
      if (sampling == "rejection")
        d_pool <- apply(genotype.pool$samples, 1, hamming.dist, y=genotype)
    }

    prob.empirical <-
      probY.empirical(N=100000, poset, lambdas, lambda.s, genotype, eps=eps)
    prob.sampling <-
      prob.importance.sampling(
        genotype, L, poset, lambdas, lambda.s, eps=eps,
        sampling.time=sampling.time, sampling=sampling,
        perturb.prob=perturb.prob, version=version, dist.pool=d_pool,
        Tdiff.pool=Tdiff_pool, seed=seeds[i])
    return(c(prob.empirical, prob.sampling))
  }

  if (!is.null(outdir)) {
    outdir <- file.path(outdir, paste("L", L, "_", sampling, sep=""))
    if (!dir.exists(outdir))
      dir.create(outdir)

    outname <- file.path(outdir, paste("probability_Y_empirical_vs_sampling",
                                       outname, ".pdf", sep=""))
    cat("Saving plot at \'", outname,"\'\n", sep="")
    if (one.genotype) {
      xlab <- expression(P(Y))
      ylab <- ""
    } else {
      xlab <- expression(P[empirical](Y))
      ylab <- expression(widehat(P)(Y))
    }

    pl.empirical_vs_sampling(
      empirical=probs[, 1], sampling=probs[, 2], xlab=xlab, ylab=ylab, N=N,
      one.genotype=one.genotype, outname=outname, binwidth=binwidth)
  }

  return(list("empirical"=probs[, 1], "sampling"=probs[, 2]))
}

#' @inheritParams prob.empirical_vs_sampling
tdiff.empirical_vs_sampling <- function(
  events, L, poset, lambdas, lambda.s, eps, rep=NULL, one.genotype=FALSE,
  sampling.times=NULL, sampling=c('forward', 'add-remove', 'rejection'),
  perturb.prob=0.8, version="3", genotype.pool=NULL, outdir=NULL, outname="",
  binwidth=0.01, seed=NULL) {

  # NOTE: If 'forward' sampling is employed, sampling times are not used.
  sampling <- match.arg(sampling)

  if (!is.null(seed))
    set.seed(seed)

  if (one.genotype) {
    if (length(sampling.times) > 1) {
      warning("Only one sampling time was expected. First entry of vector ",
              "\'sampling.times\' is used.")
      sampling.times <- sampling.times[1]
    }
    genotype <- events
    sampling.time <- sampling.times
    p <- length(genotype)
    N <- rep
  } else {
    N <- nrow(events)
    p <- ncol(events)
  }

  time.diff.empirical <- matrix(0, ncol=p, nrow=N)
  time.diff.sampling  <- matrix(0, ncol=p, nrow=N)

  for (i in 1:N) {
    if (!one.genotype) {
      genotype <- events[i, ]
      sampling.time <- sampling.times[i]
    }
    time.diff.empirical[i, ] <-
      tdiff.empirical(N=100000, poset, lambdas, lambda.s, genotype, eps=eps)
    time.diff.sampling[i, ] <-
      tdiff.importance.sampling(genotype, L=L, poset, lambdas, lambda.s, eps,
                                sampling.time=sampling.time, sampling=sampling,
                                perturb.prob=perturb.prob, version=version)
  }

  if (!is.null(outdir)) {
    outdir <- file.path(outdir, paste("L", L, "_", sampling, sep=""))
    if (!dir.exists(outdir))
      dir.create(outdir)

    for (j in 1:p) {
      pl.name <- file.path(outdir, paste("time_diff_empirical_vs_sampling",
                                         outname, "_j", j, ".pdf", sep=""))

      if (one.genotype) {
        xlab = expression(Z[j])
        ylab = ""
      } else {
        xlab = substitute(paste(Z[j]), list(j=j))
        ylab = substitute(paste(widehat(Z)[j]), list(j=j))
      }

      pl.empirical_vs_sampling(empirical=time_diff_empirical[, j],
                               sampling=time.diff.sampling[, j],
                               xlab=xlab, ylab=ylab, N=N,
                               one.genotype=one.genotype, outname=pl.name,
                               binwidth=binwidth)
    }
  }

  return(list("empirical"=time.diff.empirical, "sampling"=time.diff.sampling))
}

#' @inheritParams prob.empirical_vs_sampling
dist.empirical_vs_sampling <- function(
  events, L, poset, lambdas, lambda.s, eps, rep=NULL, one.genotype=FALSE,
  sampling.times=NULL, sampling=c('forward', 'add-remove', 'rejection'),
  perturb.prob=0.8, version="3", genotype_pool=NULL, outdir=NULL, outname="",
  binwidth=0.01, seed=NULL) {

  # NOTE: If 'forward' sampling is employed, sampling times are not used.
  sampling <- match.arg(sampling)

  if (!is.null(seed))
    set.seed(seed)

  if (one.genotype) {
    if (length(sampling.times) > 1) {
      warning("Only one sampling time was expected. First entry of vector ",
              "\'sampling.times\' is used.")
      sampling.times <- sampling.times[1]
    }
    genotype <- events
    sampling.time <- sampling.times
    N <- rep
  } else {
    N <- nrow(events)
  }

  dist.empirical <- numeric(N)
  dist.sampling  <- numeric(N)

  for (i in 1:N) {
    if (!one.genotype) {
      genotype <- events[i, ]
      sampling.time <- sampling.times[i]
    }
    dist.empirical[i] <-
      dist.empirical(N=100000, poset, lambdas, lambda.s, genotype, eps=eps)
    dist.sampling[i] <-
      dist.importance.sampling(genotype, L=L, poset, lambdas, lambda.s, eps,
                               sampling.time=sampling.time, sampling=sampling,
                               perturb.prob=perturb.prob, version=version)
  }

  if (!is.null(outdir)) {
    outdir <- file.path(outdir, paste("L", L, "_", sampling, sep=""))
    if (!dir.exists(outdir))
      dir.create(outdir)

    outname <- file.path(outdir, paste("hamming_dist_empirical_vs_sampling",
                                       outname, ".pdf", sep=""))

    if (one.genotype) {
      xlab <- expression(d(X,Y))
      ylab <- ""
    } else {
      xlab <- expression(d[empirical](X,Y))
      ylab <- expression(widehat(d)(X,Y))
    }

    pl.empirical_vs_sampling(empirical=d.empirical, sampling=d.sampling,
                             xlab=xlab, ylab=ylab, N=N,
                             one.genotype=one.genotype, outname=outname,
                             binwidth=binwidth)
  }

  return(list("empirical"=dist.empirical, "sampling"=dist.sampling))
}

#' @param empirical,sampling numeric vectors corresponding to the x and y axes,
#' respectively
#' @param xlab a title for the x axis
#' @param ylab an optional title for the y axis
#' @param truth an optional argument indicating the exact value
#' @param N an optional argument indicating the number of repetitions. This
#' option is used if \code{one.genotype} is set to \code{TRUE}
#' @param one.genotype a boolean variable indicating whether only one genotype
#' is provided
#' @param outname an optional argument indicating the output file name
#' @param binwidth an optional argument indicating the width of the histrogram
#' bins. This option is used if \code{one.genotype} is set to \code{TRUE}
pl.empirical_vs_sampling <- function(
  empirical, sampling, xlab, ylab="", truth=NULL, N=NULL, one.genotype=FALSE,
  outname=NULL, binwidth=0.01) {

  if (one.genotype) {
    df <- data.frame(x=c(empirical, sampling),
                     method=c(rep("empirical", N), rep("sampling", N)))
    pl <- ggplot(df, aes(x = x, fill = method))
    pl <- pl + geom_histogram(binwidth=binwidth, alpha=0.5, position="identity") +
      geom_vline(xintercept=truth, colour="#BB0000") +
      labs(x=xlab, y="Frequency") +
      theme_bw() + theme(text=element_text(size=14))

    if (is.null(outname)) {
      return(pl)
    } else {
      ggsave(outname, pl, width=4, height=2.5)
    }
  } else {
    df <- data.frame(x=empirical, y=sampling)
    pl <- ggplot(df, aes(x=x, y=y))
    pl <- pl + geom_point() +
      geom_abline(intercept=0, slope=1, colour="blue") +
      labs(x=xlab, y=ylab) +
      theme_bw() + theme(text=element_text(size=14))

    if (is.null(outname)) {
      return(pl)
    } else {
      ggsave(outname, pl, width=3, height=2)
    }
  }
}

#' @description initialization of the rate parameters
#'
#' @param obs.events a matrix of observations, where each row correponds to a
#' vector indicating whether an event has been observed (\code{1}) or not
#' (\code{0})
#' @param poset a matrix containing the cover relations
#' @param lambda.s rate of the sampling process
#' @param verbose an optional argument indicating whether to output logging
#' information
initialize.lambda <- function(obs.events, poset, lambda.s, verbose=FALSE) {
  p <- nrow(poset)
  lambda <- rep(0, p)
  for (i in 1:p) {
    parents <- which(poset[, i] == 1)
    if (length(parents) == 0) {
      aux <- sum(obs.events[, i]) / nrow(obs.events)
      lambda[i] <- aux * lambda.s / (1 - aux)
    } else {
      idxs <-
        which(apply(obs.events[, parents, drop=FALSE], 1,
                    function(x) all(x == 1)))

      if (length(idxs) == 0) {
        if (verbose) {
          cat("No enough evidence for mutation", i, "\n")
        }
        lambda[i] <- 1e-7
      } else {
        aux <- sum(obs.events[idxs, i]) / length(idxs)
        lambda[i] <- aux * lambda.s / (1 - aux)
      }
    }
  }
  return(lambda)
}

#' @title Monte Carlo Expectation Maximization
#' @description parameter estimation for the hidden conjunctive Bayesian network
#' model (H-CBN) via importance sampling
#'
#' @param poset a matrix containing the cover relations
#' @param obs.events a matrix of observations, where each row correponds to a
#' vector indicating whether an event has been observed (\code{1}) or not
#' (\code{0})
#' @param sampling.times an optional vector of sampling times per observation/
#' genotype
#' @param lambda.s rate of the sampling process. Defaults to \code{1.0}
#' @param max.iter the maximum number of EM iterations. Defaults to \code{100}
#' iterations
#' @param L number of samples to be drawn from the proposal in the E-step
#' @param sampling type of sampling scheme. OPTIONS: \code{"forward"} -
#' generate occurrence times according to current rate parameters, and, from
#' them, generate the genotypes; \code{"add-remove"} - generate genotypes from
#' observed genotypes using a two-steps proposal. First, make genotypes
#' compatible with the poset by either adding or removing mutations. Second,
#' perturb this version by adding or removing a mutation while yielding a
#' compatible observation; \code{"rejection"} - generate a pool of compatible
#' genotypes according to current rate parameters and sample \code{K}
#' observations proportional to the Hamming distance.
#' @param max.lambda.val an optional upper bound on the value of the rate
#' parameters. Defaults to \code{1e6}
#' @param perturb.prob probability of perturbing a genotype. Genotypes are
#' perturbed in order to learn the error rate, \eqn{\epsilon}. A genotype is
#' perturbed, if after drawing a random number, this is smaller than
#' \code{perturb.prob}, otherwise \eqn{X_start = X}. This option is used if
#' \code{sampling} is set to \code{"add-remove"}. Defaults to \code{0.3}
#' @param version an optional argument indicating which version of the
#' \code{"add-remove"} sampling scheme to use. Defaults to \code{3}
#' @param parallel boolean variable indicating whether sampling should be
#' executed sequentially (\code{0}) or in parallel (\code{1}). Option is used
#' if \code{sampling} is set to \code{"forward"}
#' @param verbose an optional argument indicating whether to output logging
#' information
#' @param seed seed for reproducibility
MCEM.hcbn_ <- function(
  poset, obs.events, sampling.times=NULL, lambda.s=1.0, max.iter=100,
  burn_in=0.8, L=100, sampling=c('forward', 'add-remove', 'rejection'),
  max.lambda.val=1e6, perturb.prob=0.3, version="3", parallel=TRUE,
  verbose=TRUE, seed=NULL) {

  # NOTE: If 'forward' sampling is employed, sampling times are not used.
  sampling <- match.arg(sampling)

  if (!is.null(seed))
    set.seed(seed, kind="L'Ecuyer-CMRG")

  p <- ncol(poset)       # number of events/mutations
  N <- nrow(obs.events)  # number of observations/genotypes

  if (is.null(sampling.times)) {
    avg.sampling.t <- lambda.s
  } else {
    avg.sampling.t <- mean(sampling.times)
    if (sampling == 'forward') {
      warning("Forward proposal doesn't account for sampling times")
    }
  }

  # initialize lambdas
  lambdas <-
    initialize.lambda(obs.events=obs.events, poset=poset, lambda.s=lambda.s)
  # NOTE: function initialize_lambda should be modify to account for cases where
  #       sum(obs.events[parents_i == 1, i]) is zero.
  lambdas[lambdas == 0] <- 1e-7

  compatible.obs <- compatible_genotypes(obs.events, poset)
  lambdas <-
    estimate_mutation_rates(
      poset, obs.events[compatible.obs$compatible_indexes, ],
      sampling.times, ilambda=lambdas, verbose=verbose)
  lambdas <- lambdas$lambda

  # initialize epsilon
  eps <- 1 - compatible.obs$fraction
  eps <- ifelse(eps == 0, runif(1, 0, 0.0001), eps)

  avg.eps <- 0
  avg.lambdas <- numeric(p)
  avg.llhood <- 0

  record.iter <- max(as.integer(burn_in * max.iter), 1)

  if (sampling == 'add-remove') {
    poset_trans_closed <- trans_closure(poset)
    parents <- get.parents(poset_trans_closed)
    children <- get.children(poset_trans_closed)
    ### **TODO** ###  at the moment, compatibility test is performed for each genotype
    #idx_compatible = apply(obs.events, 1, is_compatible, poset=poset)
  }

  llhood <- numeric()
  iter <- 1
  while(iter <= max.iter) {

    if (iter > 2) {
      if (abs(llhood[iter - 1] - llhood[iter - 2]) < 1e-4)
        break
    }

    # E step
    # Compute for each event j, expected.Tdiff[j] = E[T_j - max_{u \in pa(j)} T_u | Y]
    # Compute expected.dist
    if (sampling == 'forward') {

      if (!parallel) {
        expected.dist = numeric(N)
        expected.Tdiff = matrix(0, nrow=N, ncol=p)
        for (i in 1:N) {
          # Draw L samples from poset with parameters 'lambdas' and 'lambda.s'
          e.step <-
            importance.weight_(genotype=obs.events[i, ], L=L, poset=poset,
                               lambdas=lambdas, lambda.s=lambda.s, eps=eps,
                               sampling.time=NULL, sampling='forward')

          # Conditional expectation of the sufficient statistic d(X, Y)
          expected.dist[i] = sum(e.step$w * e.step$dist) / sum(e.step$w)

          # Contitional expectation of the sufficient statistic Z_j
          # Z_j = t_j - max_{u \in pa(j)} t_u
          expected.Tdiff[i, ] = colSums(e.step$w * e.step$time.differences) /
            sum(e.step$w)
        }
      } else {
        ret <- foreach(i=1:N, .combine='rbind', .packages="mccbn") %dopar% {

          # Draw L samples from poset with parameters 'lambdas' and 'lambda.s'
          e.step <-
            importance.weight_(genotype=obs.events[i, ], L=L, poset=poset,
                               lambdas=lambdas, lambda.s=lambda.s, eps=eps,
                               sampling.time=NULL, sampling='forward')
          # Conditional expectation of the sufficient statistic d(X, Y)
          expected.dist = sum(e.step$w * e.step$dist) / sum(e.step$w)

          # Contitional expectation of the sufficient statistic Z_j
          # Z_j = t_j - max_{u \in pa(j)} t_u
          expected.Tdiff = colSums(e.step$w * e.step$time.differences) /
            sum(e.step$w)
          return(c(expected.dist, expected.Tdiff))
        }
        expected.dist <- ret[, 1]
        expected.Tdiff <- ret[, -1]
        colnames(expected.Tdiff) <- NULL
      }

    } else if (sampling == 'add-remove') {

      if (verbose) cat("E-step - ", iter, "\n")
      ret <- foreach(i=1:N, .combine='rbind', .packages="mccbn") %dopar% {

        # In each iteration and for each observation, draw L samples from proposal
        # two-steps proposal: make compatible and perturb
        # 1. Make genotypes compatible by adding or removing 1's
        # 2. Perturb version of previous genotypes, but ensure it remains
        #    compatible with current poset
        e.step <-
          importance.weight_(genotype=obs.events[i, ], L=L, poset=poset,
                             lambdas=lambdas, lambda.s=lambda.s, eps=eps,
                             sampling.time=sampling.times[i],
                             sampling='add-remove', perturb.prob=perturb.prob,
                             version=version)
        # Conditional expectation of the sufficient statistic d(X, Y)
        expected.dist <- sum(e.step$w * e.step$dist) / sum(e.step$w)

        # Contitional expectation of the sufficient statistic Z_j
        # Z_j = t_j - max_{u \in pa(j)} t_u
        expected.Tdiff <- colSums(e.step$w * e.step$time.differences) /
          sum(e.step$w)

        return(c(expected.dist, expected.Tdiff))
      }
      expected.dist <- ret[, 1]
      expected.Tdiff <- ret[, -1]
      colnames(expected.Tdiff) <- NULL

    } else if (sampling == 'rejection') {

      K <- min(max(2^(p + 1), 2*L), 1e8 / (8 * p))
      genotype_pool <-
        sample_genotypes(K, poset, sampling_param=lambda.s, lambdas=lambdas)

      ret <- foreach(i=1:N, .combine='rbind', .packages="mccbn") %dopar% {

        # In each iteration and for each observation, draw L samples from proposal
        # rejection sampling: sample X genotypes according to the epsilon and
        # Hamming distance to the observed genotype, Y
        d_pool <-
          apply(genotype_pool$obs_events, 1, hamming.dist, y=obs.events[i, ])
        e.step <-
          importance.weight_(genotype=obs.events[i, ], L=L, poset=poset,
                             lambdas=lambdas, lambda.s=lambda.s, eps=eps,
                             sampling.time=NULL, sampling='rejection',
                             genotype.pool=genotype_pool$obs_events,
                             dist.pool=d_pool)
        # Conditional expectation of the sufficient statistic d(X, Y)
        expected.dist <- sum(e.step$w * d_pool[e.step$idxs]) / sum(e.step$w)

        # Contitional expectation of the sufficient statistic Z_j
        # Z_j = t_j - max_{u \in pa(j)} t_u
        expected.Tdiff <-
          colSums(e.step$w * genotype_pool$T_events[e.step$idxs, ]) /
          sum(e.step$w)

        return(c(expected.dist, expected.Tdiff))
      }
      expected.dist <- ret[, 1]
      expected.Tdiff <- ret[, -1]
      colnames(expected.Tdiff) <- NULL
    }
    if (any(is.na(expected.Tdiff))) {
      cat("At iteration", iter, "unexpected value for expected time differences\n")
      save(expected.Tdiff,
           file=paste("expected.Tdiff_p", p, "_iter", iter, ".RData", sep=""))
    }
    if (any(is.na(expected.dist))) {
      cat("At iteration", iter, "unexpected value for expected distance\n")
      save(expected.dist,
           file=paste("expected.dist_p", p, "_iter", iter, ".RData", sep=""))
    }

    # M step
    eps <- mean(expected.dist) / p
    lambdas <- 1 / apply(expected.Tdiff, 2, mean)
    if (verbose) cat("M-step - ", iter, "epsilon: ", eps,  "\n")
    if (any(lambdas > max.lambda.val)) {
      idx <- which(lambdas > max.lambda.val)
      lambdas[idx] <- max.lambda.val
    }

    llhood <- c(llhood,
                complete.loglikelihood_(lambdas, eps, expected.Tdiff,
                                        expected.dist))

    if (verbose)
      cat("Iteration:", iter, "- Log-likelihood: ", llhood[iter], "\n")

    if (iter > record.iter) {
      if (verbose)
        cat("Recording parameters .. \n")
      avg.eps     <- avg.eps + eps
      avg.lambdas <- avg.lambdas + lambdas
      avg.llhood  <- avg.llhood + llhood[iter]
    }
    iter <- iter + 1
  }
  avg.eps     <- avg.eps / (max.iter - record.iter)
  avg.lambdas <- avg.lambdas / (max.iter - record.iter)
  avg.llhood  <- avg.llhood / (max.iter - record.iter)
  return(list("lambdas"=lambdas, "eps"=eps, "llhood"=llhood,
              "avg.lambdas"=avg.lambdas, "avg.eps"=avg.eps,
              "avg.llhood"=avg.llhood))
}

#' @title Observed Log-Likelihood
#' @description compute the observed log-likelihood
#'
#' @param obs a matrix containing observations or genotypes, where each row
#' correponds to a genotype vector whose entries indicate whether an event has
#' been observed (\code{1}) or not (\code{0})
#' @param poset a matrix containing the cover relations
#' @param lambda a vector of the rate parameters
#' @param eps error rate
#' @param times an optional vector of sampling times per observation
#' @param L number of samples to be drawn from the proposal
#' @param sampling type of sampling scheme. OPTIONS: \code{"forward"},
#' \code{"add-remove"} or \code{"rejection"}
#' @param version an integer indicating which version of the
#' \code{"add-remove"} sampling scheme to use
#' @param genotype.pool an optional matrix containing sampled genotypes
#' @param Tdiff.pool Expected time differences for the genotype pool. This
#' option is used if \code{sampling} is set to \code{"add-remove"}
#' @param lambda.s rate of the sampling process. Defaults to \code{1.0}
#' @param thrds number of threads for parallel execution
#' @param seed seed for reproducibility
obs.loglikelihood <- function(
  obs, poset, lambda, eps, times=NULL, L,
  sampling=c('forward', 'add-remove', 'rejection'), version,
  genotype.pool=matrix(0L), Tdiff.pool=matrix(0), lambda.s=1.0, thrds=1L,
  seed=NULL) {

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
  if (!is.logical(genotype.pool))
    genotype.pool <-
      matrix(as.logical(genotype.pool), nrow=N, ncol=ncol(genotype.pool))

  if (is.null(seed))
    seed <- sample.int(3e4, 1)


  .Call("_obs_log_likelihood", PACKAGE = 'mccbn', obs, poset, lambda,
        eps, times, L, sampling, version, genotype.pool, Tdiff.pool, lambda.s,
        sampling.times.available, as.integer(thrds), as.integer(seed))
}

#' @title Sample genotypes
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
  }
  if (is.null(seed))
    seed <- sample.int(3e4, 1)

  .Call("_sample_genotypes", PACKAGE = 'mccbn', N, poset, lambda,
        matrix(0, nrow=N, ncol=p), sampling.times, lambda.s,
        sampling.times.available, as.integer(seed))
}
