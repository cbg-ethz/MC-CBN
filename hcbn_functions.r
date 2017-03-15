#' @description For binary genotypes, summing |x-y|
hamming_dist <- function(x, y) {
  dist = sum(abs(x - y))
  #dist = sum(x != y)
  return(as.integer(dist))
}

possible_genotypes <- function(p) {
  m = 2^p
  genotypes = matrix(0, nrow=m, ncol=p)
  for (i in 1:m) {
    genotypes[i, ] = gen_genotype(p, m - i + 1)
  }
  return(genotypes)
}

gen_genotype <- function(p, idx) {
  g = numeric(p)
  for(i in p:1) {
    g[i] = idx %% 2
    idx = idx %/% 2
  }
  return(g)
}

#' @description write poset to file in the format expected by H-CBN
write_poset <- function(poset, filename, datadir) {

  p <- nrow(poset)
  outfile = file.path(datadir, paste(filename, ".poset", sep=""))
  #  write header 
  write(p, outfile)
  
  relations = which(poset == 1)
  i = relations %% p
  j = (relations %/% p) + 1
  write.table(cbind(i, j), outfile,row.names=FALSE, col.names=FALSE,
              append=TRUE)
  
  #  write footer 
  write(0, outfile, append=TRUE)
}

#' @description Compute the complete-data log-likelihood or, equivalently, the
#' hidden log-likelihood
complete_log_likelihood <- function(lambdas, Tdiff, dist, eps) {

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

obs_log_likelihood <- function(
  obs_events, poset, lambdas, lambda_s, eps, sampling_times=NULL, L=1000,
  sampling=c('naive', 'add-remove', 'backward'), perturb_prob=0.8, version="3",
  genotype_pool=NULL, exact=FALSE, seed=NULL) {
  
  sampling = match.arg(sampling)
  
  if (!is.null(seed)) {
    set.seed(seed, kind="L'Ecuyer-CMRG")
  }
  
  # p: number of events / mutations
  # N: number of observations / genotypes
  if (exact) {
    ### **TODO** ### this is not correct
    p = length(lambdas)
    N = nrow(obs_events)
    llhood = 0 
    for (j in 1:p) {
      mutated = sum(obs_events[, j] == 1)
      Pr_X_0 = lambda_s/(lambda_s + lambdas[j])
      Pr_X_1 = lambdas[j]/(lambda_s + lambdas[j])
      llhood = llhood + mutated * log(eps * Pr_X_0 + (1-eps) * Pr_X_1)
      llhood = llhood + (N - mutated) * log((1-eps) * Pr_X_0 + eps * Pr_X_1)
    } 
  } else {
    prob_Y = foreach(i=1:N, .combine='c', .packages="mccbn") %dopar% {
      if (sampling == "backward") {
        d_pool = apply(genotype_pool, 1, hamming_dist, y=obs_events[i, ])
      }
      prob =
        prob_importance_sampling(
          genotype=obs_events[i, ], L=L, poset=poset, lambdas=lambdas,
          lambda_s=lambda_s, eps=eps, sampling_time=sampling_times[i],
          sampling=sampling, perturb_prob=perturb_prob, version=version,
          d_pool=d_pool)
      return(prob)
    }
    
    llhood = sum(log(prob_Y))
  }
  
  return(llhood)
}

#' @description Compute probability of genotype according to the model
#' empirically
probY_empirical <- function(N, poset, lambdas, lambda_s, genotype, eps) {
  simGenotypes <-
    sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas,
                     eps=eps)
  return(sum(apply(simGenotypes$obs_events, 1, 
                   function(x, y) all(x == y), y=genotype)) / N)
}

#' @description Compute the occurrence times expirically for a given genotype
tdiff_empirical <- function(N, poset, lambdas, lambda_s, genotype, eps) {

  simGenotypes <-
    sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas,
                     eps=eps)
  idx <- which(apply(simGenotypes$obs_events, 1, function(x, y) all(x == y),
                    y=genotype))
  
  return(apply(simGenotypes$T_events[idx, ], 2, mean))
}

#' @description Compute the average distance between genotype Y (subject to
#' noise) and N possible true genotypes, X
dist_empirical <- function(N, poset, lambdas, lambda_s, genotype, eps) {

  p <- ncol(poset)
  simGenotypes <-
    sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas,
                     eps=eps)
  idx <-
    which(apply(simGenotypes$obs_events, 1, function(x, y) all(x == y),
                y=genotype))

  return(sum(
    apply(simGenotypes$hidden_genotypes[idx, ], 1, hamming_dist, y=genotype)) /
      length(idx))
}

get_parents <- function(poset) {
  parents <- apply(poset, 2, function(x) which(x == 1))
  # Following step is need for the empty poset
  if (length(parents) == 0) {
    parents = vector("list", length=nrow(poset))
  }
  return(parents)
}

get_childreen <- function(poset) {
  childreen <- apply(poset, 1, function(x) which(x == 1))
  # Following step is need for the empty poset
  if (length(childreen) == 0) {
    childreen = vector("list", length=nrow(poset))
  }
  return(childreen)
}

#' @description Add events - look for parent nodes which have not occurred, but
#' any of their childreen had, and insert a 1 (i.e., replace 0 in parental node
#' by 1)
#'
#' @param genotype binary vector indicating whether an event has been observed
#' (1) or not (0)
#' @param childreen computed on the transitively closed poset
add_relations <- function(genotype, childreen) {
  
  p = length(genotype)
  return(sapply(1:p, 
                function(i, obs, c) 
                  ifelse(obs[i] == 0, ifelse(any(obs[c[[i]]] == 1), 1, obs[i]),
                         obs[i]), obs=genotype, c=childreen))
}

#' @description Remove events - look for child nodes which have occurred, but
#' if any of their parents have not occurred remove them (i.e., replace 1 in
#' child node by 0)
#'
#' @param genotype binary vector indicating whether an event has been observed
#' (1) or not (0)
#' @param parents computed on the transitively closed poset
remove_relations <- function(genotype, parents) {
  
  p = length(genotype)
  return(sapply(1:p, 
                function(i, obs, p) 
                  ifelse(obs[i] == 1, 
                         ifelse(any(obs[p[[i]]] == 0), 0, obs[i]), obs[i]),
                obs=genotype, p=parents))
}

get_possible_moves <- function(genotype, parents, childreen) {
  ### **TODO** ### vectorize get_possible_moves
  # Get the size of the Exit set - mutations that can happen next
  set_size_add = 0
  idxs_add = NA
  idx = which(genotype == 0)
  if (length(idx) > 0) {
    idxs_add = idx[sapply(idx, function(i, x, p) all(x[p[[i]]] == 1),
                          x=genotype, p=parents)]
    set_size_add = length(idxs_add)
  }

  # Get the size of the set of mutations that happened last - can be removed
  # and the genotype remains compatible with the poset
  set_size_remove = 0
  idxs_remove = NA
  idx = which(genotype == 1)
  if (length(idx) > 0) {
    idxs_remove = idx[sapply(idx, function(i, x, c) all(x[c[[i]]] == 0), 
                             x=genotype, c=childreen)]
    set_size_remove = length(idxs_remove)
  }
  
  return(list("set_size"= c("add"=set_size_add, "remove"=set_size_remove),
              "idxs_add"=idxs_add, "idxs_remove"=idxs_remove))
}

#' @description Perturb genotypes. Draw (uniformily) a random number between 0
#' and 1. If this number is < than perturb_prob, then perturb the genotype by
#' adding or removing an observation, but ensuring observation remains
#' compatible
#'
#' @param events matrix of genotypes, each row correponds to a vector
#' indicating whether an event has been observed (1) or not (0)
#' @param parents computed on the transitively closed poset
#' @param childreen computed on the transitively closed poset
#' @return List with folliwing named arguments: (1) X - perturbed version of
#' compatible genotypes, remains compatible with the current poset,
#' (2) option_set_size size of the exit set or set of observations which can be
#' removed. NA is returned when the sample was not perturbed
perturb <- function(events, parents, childreen, perturb_prob=0.8, 
                    compatible=FALSE, version=c("1", "2", "3")) {
  
  version = match.arg(version)
  
  L = nrow(events)  # number of samples
  p = ncol(events)  # number of events/mutations
  new_obs = events
  option_set_size = rep(NA, L) 
  
  perturb = ifelse(runif(L, 0, 1) < perturb_prob, 1, 0)
  idxs_perturb = which(perturb == 1)
  
  # When observation Y was already compatible, all rows are equivalent
  if (compatible || var(as.vector(events)) == 0) {
    # Get indices of mutations that can be either added or removed
    idxs = get_possible_moves(events[1, ], parents, childreen)
  }
  # Deal with special cases: wild type and resistant type
  if (all(events == 0)) {
    if (version == "1") {
      add = runif(length(idxs_perturb), 0, 1) < 0.5
    } else if (version == "2" || version == "3") {
      add = rep(TRUE, length(idxs_perturb))
    }
    mutation_idx = sample(idxs$idxs_add, sum(add), replace=TRUE)
    option_set_size[idxs_perturb[add]] = idxs$set_size["add"]
    if (sum(add) == 1) {
      new_obs[idxs_perturb[add], mutation_idx] = 1
    } else {
      ### **TODO** ### vectorize it
      for (i in 1:sum(add)) {
        new_obs[idxs_perturb[add][i], mutation_idx[i]] = 1
      }
    }
  } else if (all(events == 1)) {
    if (version == "1") {
      add = runif(length(idxs_perturb), 0, 1) < 0.5
    } else if (version == "2" || version == "3") {
      add = rep(FALSE, length(idxs_perturb))
    }
    mutation_idx = sample(idxs$idxs_remove, sum(!add), replace=TRUE)
    option_set_size[idxs_perturb[!add]] = idxs$set_size["remove"]
    if (sum(!add) == 1) {
      new_obs[idxs_perturb[!add], mutation_idx] = 0
    } else {
      ### **TODO** ### vectorize it
      for (i in 1:sum(!add)) {
        new_obs[idxs_perturb[!add][i], mutation_idx[i]] = 0
      }
    }
  } else {
    for (i in idxs_perturb) {

      # Get indices of mutations that can be either added or removed
      if (!compatible) {
        idxs = get_possible_moves(events[i, ], parents, childreen)
      }

      if (version == "1") {
        # Add or remove one observation, both moves have equal probability.
        # Draw a random number between 0 and 1. If this number is > 0.5, then
        # add an observation
        add = runif(1, 0, 1) < 0.5
        if (add) {
          mutation_idx <-
            ifelse(idxs$set_size["add"] > 1, sample(idxs$idxs_add, 1), idxs$idxs_add)
          option_set_size[i] = idxs$set_size["add"]
          new_obs[i, mutation_idx] = 1
        } else {
          mutation_idx <-
            ifelse(idxs$set_size["remove"] > 1, sample(idxs$idxs_remove, 1),
                   idxs$idxs_remove)
          option_set_size[i] = idxs$set_size["remove"]
          new_obs[i, mutation_idx] = 0
        }
      } else if (version == "2" || version == "3") {
        option_set_size[i] = sum(idxs$set_size)

        # Choose one index randomly
        idx = sample.int(option_set_size[i], 1)
        mutation_idx = na.omit(c(idxs$idxs_add, idxs$idxs_remove))[idx]

        # Check whether move corresponds to an adding move
        add = ifelse(idx <= idxs$set_size["add"], TRUE, FALSE)

        new_obs[i, mutation_idx] = ifelse(add, 1, 0)
      }
    }
  }
  return(list("X"=new_obs, "option_set_size"=option_set_size))
}

#' @description Draw L samples from the proposal. Two-steps proposal: (1) make
#' compatible and (2) perturb genotype
#'
#' @param genotype binary vector indicating whether an event has been observed
#' (1) or not (0)
#' @param L number of samples
#' @param parents computed on the transitively closed poset
#' @param childreen computed on the transitively closed poset
#' @param compatible variable indicating whether the genotype is compatible (1)
#' with current poset or not (0)
#' @param perturb_prob genotypes are perturbed in order to learn epsilon. A
#' genotype is perturbed, if after drawing a random number, this is less than
#' perturb_prob, otherwise X_start = X
draw_samples <- function(genotype, L, parents, childreen, compatible, eps,
                         perturb_prob=0.8, version=c("1", "2", "3")) {

  version = match.arg(version)
  
  p = length(genotype)  # number of events/mutations
  
  if (compatible) {
    # if genotype is compatible with poset, we can inmediately go to step 2:
    # perturb
    compatible_obs = matrix(genotype, nrow=L, ncol=p, byrow=TRUE)
    make_compatible_prob = NA
  } else {
    # if genotype is NOT compatible with poset, we need to generate L compatible
    # versions by adding or removing observations
    compatible_by_adding = add_relations(genotype, childreen)
    compatible_by_removing = remove_relations(genotype, parents)
    
    if (version == "1" | version == "2") {
      # Both moves, add and remove, have equal probability
      # Move add (remove) wouldn't be applicable if genotype is the resistant
      # type (wild type). However, such genotypes are always compatible with
      # the poset
      add = ifelse(runif(L, 0, 1) < 0.5, 1, 0)
      make_compatible_prob = rep(0.5, L)
    } else if (version == "3") {
      # Choose to add or remove observations according to the number of
      # modifications
      dist_add = hamming_dist(genotype, compatible_by_adding)
      dist_remove = hamming_dist(genotype, compatible_by_removing)
      add_prob = eps^dist_add / (eps^dist_add + eps^dist_remove)
      add = ifelse(runif(L, 0, 1) < add_prob, 1, 0)
      make_compatible_prob = ifelse(add == 1, add_prob, 1 - add_prob)
    }
    compatible_obs = matrix(0, L, p)
    
    if (all(add == 1)) {
      compatible_obs <- matrix(compatible_by_adding, nrow=L, ncol=p, byrow=TRUE)
    } else if (all(add == 0)) {
      compatible_obs <-
        matrix(compatible_by_removing, nrow=L, ncol=p, byrow=TRUE)
    } else {
      idxs_add = which(add == 1)
      compatible_obs[idxs_add, ] = matrix(compatible_by_adding, nrow=sum(add), 
                                          ncol=p, byrow=TRUE)
      compatible_obs[-idxs_add, ] = matrix(compatible_by_removing, nrow=L-sum(add), 
                                           ncol=p, byrow=TRUE)
    }
  }
  samples <-
    perturb(compatible_obs, parents, childreen, perturb_prob, compatible,
            version=version)
  return(list("samples"=samples, "make_compatible_prob"=make_compatible_prob))
}

rtexp <- function(x, rate) {
  rand_num = runif(1, 0, 1)
  Z = ifelse(x > 0, -log(1 - rand_num * (1 - exp(-rate * x))) / rate, 0) 
  return(Z)
}

sample_mutation_times <- function(genotype, poset, lambdas, lambda_s=1.0,
                                  sampling_time=NULL) {
  
  p = length(genotype)  # number of events/mutations
  parents = get_parents(poset)
  if (is.null(sampling_time)) {
    sampling_time = rexp(1, rate=lambda_s)
  }
  Tdiff = numeric(p)
  Tsum = numeric(p)
  dens = 0
  topo_path = my.topological.sort(poset)
  
  for (j in topo_path) {
    max_parents_t = 0
    if (length(parents[[j]]) > 0) {
      for (u in 1:length(parents[[j]])) {
        if(Tsum[parents[[j]][u]] > max_parents_t) {
          max_parents_t = Tsum[parents[[j]][u]]
        }
      }
    }
    
    if (genotype[j] == 1) {
      # If mutation is observed, Z ~ TExp(lambda, 0, sampling_time - max_parents_t)
      # When a lot of mutations have occurred max_parents_t approaches the 
      # sampling_time. Due to numerical precision problems, 
      # sampling_time < max_parents_t. Therefore, we have to make sure that when it 
      # happens rtexp returns 0
      # dexp: P(Z)
      # pexp: P(Z <= sampling_time - max_parents_t)
      Z = rtexp(sampling_time - max_parents_t, rate=lambdas[j])
      Tsum[j] =  max_parents_t + Z
      dens = dens +  dexp(Z, rate=lambdas[j], log=TRUE) - 
        pexp(sampling_time - max_parents_t, rate=lambdas[j], log.p=TRUE)
    } else {
      # If mutation is not observed, Z ~ Exp(lambda)
      Z = rexp(1, rate=lambdas[j])
      Tsum[j] = max(sampling_time, max_parents_t) + Z
      dens =  dens + dexp(Z, rate=lambdas[j], log=TRUE)
    }
    Tdiff[j] = Tsum[j] - max_parents_t
  }
  return(c("density"=dens, "Tdiff"=Tdiff))
}

log_cbn_density <- function(Tdiff, rate) {
  dens = sum(log(rate)) - sum(rate * Tdiff)
  return(dens)
}


importance_weight <- function(
  genotype, L, poset, lambdas, lambda.s, eps, sampling.time=NULL,
  sampling=c('naive', 'add-remove', 'backward'), perturb.prob=0.8, version="3",
  dist.pool=NULL, seed=NULL) {
  
  sampling <- match.arg(sampling)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  p <- ncol(poset) # number of events/mutations
  
  if (sampling == 'naive') {
    
    # Generate L samples from poset with parameters 'lambdas' and 'lambda_s'
    # In particular, epsilon is zero (default value) - because the idea is to
    # generate samples of X (underlying true)
    ### **TODO** ### what to do when sampling times available. At the moment,
    # not considered to generate samples - incorporate sampling times in
    # function sample_genotypes
    if (!is.null(sampling.time)) {
      warning("Naive proposal doesn't account for sampling times")
    }
    samples <- sample_genotypes(L, poset, sampling_param=lambda.s, lambdas=lambdas)
    d <- apply(samples$hidden_genotypes, 1, hamming_dist, y=genotype)
    prob_Y_X <- eps^d * (1-eps)^(p-d)
    
    return(list("w"=prob_Y_X, "time_differences"=samples$T_events, "dist"=d))
    
  } else if (sampling == 'add-remove') {
    
    poset_trans_closed = trans_closure(poset)
    parents = get_parents(poset_trans_closed)
    childreen = get_childreen(poset_trans_closed)
    compatible = is_compatible(genotype, poset)
    # Generate L samples according to proposal
    return_list <-
      draw_samples(genotype, L, parents, childreen, compatible, eps,
                   perturb.prob, version)
    make_compatible_prob <- return_list$make_compatible_prob
    samples <- return_list$samples
    
    # Generate mutation times Tj from sample i
    Tdiff = apply(samples$X, 1, sample_mutation_times, poset=poset,
                  lambdas=lambdas, lambda_s=lambda.s,
                  sampling_time=sampling.time)
    
    log_proposal_X = Tdiff["density", ]
    Tdiff = t(tail(Tdiff, n=-1))
    
    # Hamming distance bewteen samples and genotype
    dist = apply(samples$X, 1, hamming_dist, y=genotype)
    
    # Computing log(Pr(Y|X))
    if (eps == 0) {
      # NOTE: If all observations are compatible with poset and pertub_prob == 0
      #       (which can happen because noisy observations can be compatible),
      #       then eps can be 0 and dist 0
      log_Pr_Y_X <-
        ifelse(dist == 0, 0,
               dist * log(eps + .Machine$double.eps) +
                 (p - dist) * log(1 - (eps + .Machine$double.eps)))
    } else {
      log_Pr_Y_X <- dist*log(eps) + (p - dist) * log(1 - eps)
    }
    
    # Computing log(Pr(X))
    log_Pr_X = apply(Tdiff, 1, log_cbn_density, rate=lambdas)
    
    # Computing density of the proposal - for correction
    # proposal : weight accounting for making genotype compatible + weight 
    #            accounting for choosing a possible mutation for perturbing
    #            the sample + weight accounting for sampling times from 
    #            a truncated exponential
    num_options = samples$option_set_size
    if (version == "1") {
      log_proposal_Y_X <-
        ifelse(is.na(num_options) | num_options == 0, log(1 - perturb.prob),
               log(perturb.prob) + log(0.5) + log(1/num_options))
    } else if (version == "2" || version == "3") {
      stopifnot(all(na.omit(num_options) != 0))
      log_proposal_Y_X <- ifelse(is.na(num_options), log(1 - perturb.prob),
                                log(perturb.prob) + log(1 / num_options))
    }

    if (!compatible) {
      log_proposal_Y_X = log_proposal_Y_X + log(make_compatible_prob) 
    } 
    
    log_proposal = log_proposal_Y_X + log_proposal_X
    
    importance_weight = exp(log_Pr_X + log_Pr_Y_X - log_proposal)
    return(list("w"=importance_weight, "time_differences"=Tdiff, "dist"=dist))

  } else if (sampling == 'backward') {
    ### **TODO** ### what to do when sampling times available. At the moment,
    # not considered to generate samples - incorporate sampling times in
    # function sample_genotypes
    if (is.null(dist.pool)) {
      stop("Vector of distances expected for backward sampling")
    }
    # if (is.null(genotype_pool)) {
    #   stop("Pool of genotypes are expected for backward sampling")
    # }
    K <- length(dist.pool)
    # compatible = is_compatible(genotype, poset)
    aux <- eps^(dist.pool) * (1-eps) ^ (p - dist.pool)
    # if (compatible) {
    #   genotype_pool = rbind(genotype_pool, genotype)
    #   q_prob = c(aux, 1 - sum(aux))
    #   idxs = sample(seq(1, K + 1), L, replace=TRUE, prob=q_prob)
    # } else {
    #   q_prob = aux/sum(aux)
    #   idxs = sample(seq(1, K), L, replace=TRUE, prob=q_prob)
    # }
    q_prob <- aux / sum(aux)
    idxs <- sample(seq(1, K), L, replace=TRUE, prob=q_prob)
    # Generate mutation times Tj from sample i
    # Tdiff <-
    #   apply(genotype_pool[idxs, ], 1, sample_mutation_times, poset=poset,
    #         lambdas=lambdas, lambda_s=lambda.s, sampling_time=sampling.time)
    #
    # log_proposal_X = Tdiff["density", ]
    # Tdiff = t(tail(Tdiff, n=-1))

    # Hamming distance bewteen samples and genotype
    dist <- dist.pool[idxs]

    # Computing log(Pr(Y|X))
    if (eps == 0) {
      # NOTE: If all obs ervations are compatible with poset and pertub_prob == 0
      #       (which can happen because noisy observations can be compatible),
      #       then eps can be 0 and dist 0
      log_Pr_Y_X <-
        ifelse(dist == 0, 0,
               dist * log(eps + .Machine$double.eps) +
                 (p - dist) * log(1 - (eps + .Machine$double.eps)))
    } else {
      log_Pr_Y_X <- dist * log(eps) + (p - dist) * log(1 - eps)
    }

    # Computing log(Pr(X))
    # log_Pr_X = apply(Tdiff, 1, log_cbn_density, rate=lambdas)

    # Computing density of the proposal - for correction
    log_proposal <- log(q_prob[idxs]) + log(K) #+ log_proposal_X

    # importance_weight = exp(log_Pr_X + log_Pr_Y_X - log_proposal)
    importance_weight <- exp(log_Pr_Y_X - log_proposal)
    # return(list("w"=importance_weight, "time_differences"=Tdiff, "dist"=dist))
    return(list("w"=importance_weight, "idxs"=idxs))
  }
}

#' @description Compute Pr(Y) using Monte Carlo sampling, where Y corresponds
#' to the observed genotype
#'
#' @param L number of samples
#' @param sampling type of sampling. Options: 'naive' - generating occurrence
#' times according to lambdas and from them genotypes; 'add-remove' -
#' generating genotypes from observed genotypes using a two-steps proposal.
#' First, making genotypes compatible with the poset by either adding or
#' removing observations. Second, perturbing this version by adding or removing
#' observations
prob_importance_sampling <- function(
  genotype, L, poset, lambdas, lambda_s, eps, sampling_time=NULL,
  sampling=c('naive', 'add-remove', 'backward'), perturb_prob=0.8, version="3",
  d_pool=NULL, seed=NULL) {
  # NOTE: If 'naive' sampling is employed, sampling times are not used. 
  
  sampling = match.arg(sampling)
  if (!is.null(seed)) {
    set.seed(seed, kind="L'Ecuyer-CMRG")
  }
  probs <-
    importance_weight(genotype, L, poset, lambdas, lambda_s, eps,
                      sampling.time=sampling_time, sampling=sampling,
                      perturb.prob=perturb_prob, version=version,
                      dist.pool=d_pool)
  probs <- probs$w

  return(sum(probs) / L)
}

#' @description Compute expected time differences for a given genotype using
#' Monte Carlo sampling
tdiff_importance_sampling <-
  function(genotype, L, poset, lambdas, lambda_s, eps, sampling_time=NULL,
           sampling=c('naive', 'add-remove', 'backward'), perturb_prob=0.8,
           version="3", genotype_pool=NULL, d_pool=NULL) {

  sampling = match.arg(sampling)
  importance_weights = 
    importance_weight(genotype, L, poset, lambdas, lambda_s, eps, 
                      sampling.time=sampling_time, sampling=sampling,
                      perturb.prob=perturb_prob, version=version,
                      dist.pool=d_pool)
  return(colSums(importance_weights$w * importance_weights$time_differences) / 
           sum(importance_weights$w))

}

#' @description Compute expected hamming distance between a given observations
#' (genotype) and the underlying/true genotype using Monte Carlo sampling
dist_importance_sampling <-
  function(genotype, L, poset, lambdas, lambda_s, eps, sampling_time=NULL,
           sampling=c('naive', 'add-remove', 'backward'), perturb_prob=0.8,
           version=version, genotype_pool=NULL, d_pool=NULL) {

  sampling = match.arg(sampling)
  importance_weights = 
    importance_weight(genotype, L, poset, lambdas, lambda_s, eps, 
                      sampling.time=sampling_time, sampling=sampling,
                      perturb.prob=perturb_prob, version=version,
                      dist.pool=d_pool)
  return(sum(importance_weights$dist * importance_weights$w) / 
           sum(importance_weights$w))
}

#' @param events matrix of observations or single genotype, if 'one_genotype' is TRUE
#' @param L number of samples
#' @param poset
#' @param lambda vector of exponential rates
#' @param eps
#' @param rep if only one genotype is provided, number of repetitions
#' @param one_genotype boolean variable indicating if only one genotype is provided
#' @param sampling_times vector of sampling times per event
#' @param sampling type of sampling scheme. Options: 'naive', 'add-remove' or 'backward'
#' @param perturb_prob option used if sampling is set to 'add-remove'
#' @return  returns importance weights and sufficient statistics
prob_empirical_vs_sampling <- function(
  events, L, poset, lambdas, lambda_s, eps, rep=NULL, one_genotype=FALSE,
  sampling_times=NULL, sampling=c('naive', 'add-remove', 'backward'),
  perturb_prob=0.8, version="3", genotype_pool=NULL, outdir=NULL, outname="",
  binwidth=0.01, seed=NULL) {
  # NOTE: If 'naive' sampling is employed, sampling times are not used. 
  
  sampling <- match.arg(sampling)
  
  if (!is.null(seed)) {
    set.seed(seed, kind="L'Ecuyer-CMRG")
  }
  if (one_genotype) {
    if (length(sampling_times) > 1) {
      warning("Only one sampling time was expected. First entry of vector ", 
               "\'sampling_times\' is used.")
      sampling_times <- sampling_times[1]
    }
    genotype <- events
    sampling_time <- sampling_times
    N = rep
  } else {
    N = nrow(events)
  }
  
  prob_empirical = numeric(N)
  prob_sampling  = numeric(N)
  
  probs = foreach(i=1:N, .combine='rbind', .packages="mccbn") %dopar% {
    if (!one_genotype) {
      genotype <- events[i, ]
      sampling_time <- sampling_times[i]
    }
    if (sampling == "backward") {
      d_pool <- apply(genotype_pool, 1, hamming_dist, y=genotype)
    }
    prob_empirical <-
      probY_empirical(N=100000, poset, lambdas, lambda_s, genotype, eps=eps)
    prob_sampling <-
      prob_importance_sampling(
        genotype, L, poset, lambdas, lambda_s, eps=eps,
        sampling_time=sampling_time, sampling=sampling,
        perturb_prob=perturb_prob, version=version, d_pool=d_pool)
    return(c(prob_empirical, prob_sampling))
  }
  
  if (!is.null(outdir)) {
    outdir <- file.path(outdir, paste("L", L, "_", sampling, sep=""))
    if (!dir.exists(outdir)) {
      dir.create(outdir)
    }
    outname <- file.path(outdir, paste("probability_Y_empirical_vs_sampling",
                                       outname, ".pdf", sep=""))
    cat("Saving plot at \'", outname,"\'\n", sep="")
    if (one_genotype) { 
      xlab <- expression(P(Y))
      ylab <- ""
    } else {
      xlab <- expression(P[empirical](Y))
      ylab <- expression(widehat(P)(Y))
    }
    
    pl_empirical_vs_sampling(
      empirical=probs[, 1], sampling=probs[, 2], xlab=xlab, ylab=ylab, N=N,
      one_genotype=one_genotype, outname=outname, binwidth=binwidth)
  }

  return(list("empirical"=probs[, 1], "sampling"=probs[, 2]))
}

#' @param events matrix of observations or single genotype, if 'one_genotype'
#' is TRUE
#' @param L number of samples
#' @param rep if only one genotype is provided, number of repetitions
#' @param one_genotype boolean variable indicating if only one genotype is
#' provided
#' @param sampling_times vector of sampling times per event
#' @param sampling type of sampling. OPTIONS: 'naive' or 'add-remove'
tdiff_empirical_vs_sampling <- function(
  events, L, poset, lambdas, lambda_s, eps, rep=NULL, one_genotype=FALSE,
  sampling_times=NULL, sampling=c('naive', 'add-remove', 'backward'),
  perturb_prob=0.8, version="3", genotype_pool=NULL, outdir=NULL, outname="",
  binwidth=0.01, seed=NULL) {
  # NOTE: If 'naive' sampling is employed, sampling times are not used. 
  
  sampling = match.arg(sampling)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (one_genotype) {
    if (length(sampling_times) > 1) {
      warning("Only one sampling time was expected. First entry of vector ", 
              "\'sampling_times\' is used.")
      sampling_times = sampling_times[1]
    }
    genotype = events
    sampling_time = sampling_times
    p = length(genotype)
    N = rep
  } else {
    N = nrow(events)
    p = ncol(events)
  }
  
  time_diff_empirical = matrix(0, ncol=p, nrow=N)
  time_diff_sampling  = matrix(0, ncol=p, nrow=N)
  
  for (i in 1:N) {
    if (!one_genotype) {
      genotype = simulated_obs$obs_events[i, ]
      sampling_time = sampling_times[i]
    }
    time_diff_empirical[i, ] <-
      tdiff_empirical(N=100000, poset, lambdas, lambda_s, genotype, eps=eps)
    time_diff_sampling[i, ] <-
      tdiff_importance_sampling(genotype, L=L, poset, lambdas, lambda_s, eps,
                                sampling_time=sampling_time, sampling=sampling,
                                perturb_prob=perturb_prob, version=version)
  }
  
  if (!is.null(outdir)) {
    outdir = file.path(outdir, paste("L", L, "_", sampling, sep=""))
    if (!dir.exists(outdir)) {
      dir.create(outdir)
    }
    
    for (j in 1:p) {
      pl_name <- file.path(outdir, paste("time_diff_empirical_vs_sampling",
                                         outname, "_j", j, ".pdf", sep=""))

      if (one_genotype) {
        xlab = expression(Z[j])
        ylab = ""
      } else {
        xlab = substitute(paste(Z[j]), list(j=j))
        ylab = substitute(paste(widehat(Z)[j]), list(j=j))
      }

      pl_empirical_vs_sampling(empirical=time_diff_empirical[, j],
                               sampling=time_diff_sampling[, j],
                               xlab=xlab, ylab=ylab, N=N,
                               one_genotype=one_genotype, outname=pl_name,
                               binwidth=binwidth)

    }
  }
  
  return(list("empirical" = time_diff_empirical, "sampling" = time_diff_sampling))
}

#' @param events matrix of observations or single genotype, if 'one_genotype'
#' is TRUE 
#' @param L number of samples
#' @param rep if only one genotype is provided, number of repetitions
#' @param one_genotype boolean variable indicating if only one genotype is
#' provided
#' @param sampling_times vector of sampling times per event
#' @param sampling type of sampling. OPTIONS: 'naive' or 'add-remove'.
dist_empirical_vs_sampling <- function(
  events, L, poset, lambdas, lambda_s, eps, rep=NULL, one_genotype=FALSE,
  sampling_times=NULL, sampling=c('naive', 'add-remove', 'backward'),
  perturb_prob=0.8, version="3", genotype_pool=NULL, outdir=NULL, outname="",
  binwidth=0.01, seed=NULL) {
  # NOTE: If 'naive' sampling is employed, sampling times are not used. 

  sampling = match.arg(sampling)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (one_genotype) {
    if (length(sampling_times) > 1) {
      warning("Only one sampling time was expected. First entry of vector ", 
              "\'sampling_times\' is used.")
      sampling_times = sampling_times[1]
    }
    genotype = events
    sampling_time = sampling_times
    N = rep
  } else {
    N = nrow(events)
  }
  
  d_empirical = numeric(N)
  d_sampling  = numeric(N)
  
  for (i in 1:N) {
    if (!one_genotype) {
      genotype = simulated_obs$obs_events[i, ]
      sampling_time = sampling_times[i]
    }
    d_empirical[i] <-
      dist_empirical(N=100000, poset, lambdas, lambda_s, genotype, eps=eps)
    d_sampling[i] <-
      dist_importance_sampling(genotype, L=L, poset, lambdas, lambda_s, eps,
                               sampling_time=sampling_time, sampling=sampling,
                               perturb_prob=perturb_prob, version=version)
  }
  
  if (!is.null(outdir)) {
    outdir = file.path(outdir, paste("L", L, "_", sampling, sep=""))
    if (!dir.exists(outdir)) {
      dir.create(outdir)
    }
    outname <- file.path(outdir, paste("hamming_dist_empirical_vs_sampling",
                                       outname, ".pdf", sep=""))

    if (one_genotype) {
      xlab = expression(d(X,Y))
      ylab = ""
    } else {
      xlab = expression(d[empirical](X,Y))
      ylab = expression(widehat(d)(X,Y))
    }

    pl_empirical_vs_sampling(empirical=d_empirical, sampling=d_sampling,
                             xlab=xlab, ylab=ylab, N=N,
                             one_genotype=one_genotype, outname=outname,
                             binwidth=binwidth)
  }
  
  return(list("empirical" = d_empirical, "sampling" = d_sampling))
}


pl_empirical_vs_sampling <- function(
  empirical, sampling, xlab, ylab="", truth=NULL, N=NULL, one_genotype=FALSE,
  outname=NULL, binwidth=0.01) {
  
  if (one_genotype) {
    df = data.frame(x=c(empirical, sampling), 
                    method=c(rep("empirical", N), rep("sampling", N)))
    pl = ggplot(df, aes(x = x, fill = method))
    pl = pl + geom_histogram(binwidth=binwidth, alpha=0.5, position="identity") + 
      geom_vline(xintercept=truth, colour="#BB0000") +
      labs(x=xlab, y="Frequency") + 
      theme_bw() + theme(text=element_text(size=14))
    
    if (is.null(outname)) {
      return(pl)
    } else {
      ggsave(outname, pl, width=4, height=2.5)
    }
  } else {
    df = data.frame(x=empirical, y=sampling)
    pl = ggplot(df, aes(x=x, y=y))
    pl = pl + geom_point() + 
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

initialize_lambda2 <- function(obs.events, poset, lambda.s, verbose=FALSE) {
  p <- nrow(poset)
  lambda <- rep(0, p)
  for (i in 1:p) {
    parents = which(poset[, i] == 1)
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

#' @param obs_events matrix of observations, each row correponds to a vector
#' indicating whether an event has been observed (1) or not (0)
#' @param sampling_times vector of sampling times per observation/genotype
#' @param L number of samples to be drawn from the proposal
#' @param sampling type of sampling. Options: 'naive', 'add-remove', 'backward'
#' @param perturb_prob genotypes are perturbed in order to learn epsilon. A
#' genotype is perturbed, if after drawing a random number, this is <
#' perturb_prob, otherwise X_start = X. Option is used if sampling is set to
#' 'add-remove'
#' @param parallel boolean variable indicating whether sampling should be
#' executed sequentially (0) or in parallel (1). Option is used if sampling is
#' set to 'naive'.
MCEM_hcbn <- function(
  poset, obs_events, sampling_times=NULL, lambda_s=1.0, max_iter=100,
  burn_in=0.8, L=100, sampling=c('naive', 'add-remove', 'backward'),
  max_lambda_val=1e6, perturb_prob=0.8, version="3", parallel=TRUE,
  verbose=TRUE, seed=NULL) {
  # NOTE: If 'naive' sampling is employed, sampling times are not used.
  
  sampling <- match.arg(sampling)
  
  if (!is.null(seed)) {
    set.seed(seed, kind="L'Ecuyer-CMRG")
  }
  
  p = ncol(poset)       # number of events/mutations
  N = nrow(obs_events)  # number of observations/genotypes
  
  if (is.null(sampling_times)) {
    avg_sampling_t = lambda_s
    
  } else {
    avg_sampling_t = mean(sampling_times)
    if (sampling == 'naive') {
      warning("Naive proposal doesn't account for sampling times")
    }
  }
  
  # initialize lambdas
  lambdas <-
    initialize_lambda2(obs.events=obs_events, poset=poset, lambda.s=lambda_s)
  # NOTE: function initialize_lambda should be modify to account for cases where
  #       sum(obs_events[parents_i == 1, i]) is zero. 
  lambdas[lambdas == 0] <- 1e-7

  compatible.obs <- compatible_genotypes(obs_events, poset)
  lambdas <-
    estimate_mutation_rates(
      poset, obs_events[compatible.obs$compatible_indexes, ],
      sampling_times, ilambda=lambdas, verbose=verbose)
  lambdas <- lambdas$lambda

  # initialize epsilon
  eps <- 1 - compatible.obs$fraction
  eps <- ifelse(eps == 0, runif(1, 0, 0.0001), eps)
  
  avg_eps = 0 
  avg_lambdas = numeric(p)
  avg_llhood = 0
  
  record_iter = max(as.integer(burn_in * max_iter), 1)
  
  if (sampling == 'add-remove') {
    poset_trans_closed = trans_closure(poset)
    parents = get_parents(poset_trans_closed)
    childreen = get_childreen(poset_trans_closed)
    ### **TODO** ###  at the moment, compatibility test is performed for each genotype
    #idx_compatible = apply(obs_events, 1, is_compatible, poset=poset)
  }
  
  llhood = numeric()

  iter = 1
  while(iter <= max_iter) {
    
    if (iter > 2) {
      if (abs(llhood[iter - 1] - llhood[iter - 2]) < 1e-4)
        break
    }
    # E step
    # Compute for each event j, expected_Tdiff[j] = E[T_j - max_{u \in pa(j)} T_u | Y]
    # Compute expected_dist
    
    if (sampling == 'naive') {
      
      if (!parallel) {
        expected_dist = numeric(N)
        expected_Tdiff = matrix(0, nrow=N, ncol=p)
        for (i in 1:N) {
          # Draw L samples from poset with parameters 'lambdas' and 'lambda_s'
          e_step <-
            importance_weight(genotype=obs_events[i, ], L=L, poset=poset,
                              lambdas=lambdas, lambda.s=lambda_s, eps=eps,
                              sampling.time=NULL, sampling='naive')

          # Conditional expectation of the sufficient statistic d(X, Y)
          expected_dist[i] = sum(e_step$w * e_step$dist) / sum(e_step$w)
          
          # Contitional expectation of the sufficient statistic Z_j
          # Z_j = t_j - max_{u \in pa(j)} t_u
          expected_Tdiff[i, ] = colSums(e_step$w * e_step$time_differences) /
            sum(e_step$w)
        }
      } else {
        ret = foreach(i=1:N, .combine='rbind', .packages="mccbn") %dopar% {
          
          # Draw L samples from poset with parameters 'lambdas' and 'lambda_s'
          e_step <-
            importance_weight(genotype=obs_events[i, ], L=L, poset=poset,
                              lambdas=lambdas, lambda.s=lambda_s, eps=eps,
                              sampling.time=NULL, sampling='naive')
          # Conditional expectation of the sufficient statistic d(X, Y)
          expected_dist = sum(e_step$w * e_step$dist) / sum(e_step$w)
          
          # Contitional expectation of the sufficient statistic Z_j
          # Z_j = t_j - max_{u \in pa(j)} t_u
          expected_Tdiff = colSums(e_step$w * e_step$time_differences) /
            sum(e_step$w)
          return(c(expected_dist, expected_Tdiff))
        }
        expected_dist = ret[, 1]
        expected_Tdiff = ret[, -1]
        colnames(expected_Tdiff) = NULL 
      }
      
    } else if (sampling == 'add-remove') {
      
      if (verbose) cat("E-step - ", iter, "\n")
      ret = foreach(i=1:N, .combine='rbind', .packages="mccbn") %dopar% {
        
        # In each iteration and for each observation, draw L samples from proposal
        # two-steps proposal: make compatible and perturb
        # 1. Make genotypes compatible by adding or removing 1's
        # 2. Perturbe version of previous genotypes, but ensure it remains
        #    compatible with current poset
        e_step <-
          importance_weight(
            genotype=obs_events[i, ], L=L, poset=poset, lambdas=lambdas,
            lambda.s=lambda_s, eps=eps, sampling.time=sampling_times[i],
            sampling='add-remove', perturb.prob=perturb_prob, version=version)
        # Conditional expectation of the sufficient statistic d(X, Y)
        expected_dist = sum(e_step$w * e_step$dist) / sum(e_step$w)
        
        # Contitional expectation of the sufficient statistic Z_j
        # Z_j = t_j - max_{u \in pa(j)} t_u
        expected_Tdiff = colSums(e_step$w * e_step$time_differences) /
          sum(e_step$w)

        return(c(expected_dist, expected_Tdiff))

      }
      expected_dist = ret[, 1]
      expected_Tdiff = ret[, -1]
      colnames(expected_Tdiff) = NULL

    } else if (sampling == 'backward') {

      K = min(max(2^(p + 1), 2*L), 1e8 / (8 * p))
      genotype_pool <-
        sample_genotypes(K, poset, sampling_param=lambda_s, lambdas=lambdas)

      ret = foreach(i=1:N, .combine='rbind', .packages="mccbn") %dopar% {

        # In each iteration and for each observation, draw L samples from proposal
        # backward proposal: sample X genotypes according to the epsilon and
        # Hamming distance to the observed genotype, Y
        d_pool <-
          apply(genotype_pool$obs_events, 1, hamming_dist, y=obs_events[i, ])
        e_step <-
          importance_weight(genotype=obs_events[i, ], L=L, poset=poset,
                            lambdas=lambdas, lambda.s=lambda_s, eps=eps,
                            sampling.time=NULL, sampling='backward',
                            genotype.pool=genotype_pool$obs_events,
                            dist.pool=d_pool)
        # Conditional expectation of the sufficient statistic d(X, Y)
        expected_dist <- sum(e_step$w * d_pool[e_step$idxs]) / sum(e_step$w)

        # Contitional expectation of the sufficient statistic Z_j
        # Z_j = t_j - max_{u \in pa(j)} t_u
        expected_Tdiff <-
          colSums(e_step$w * genotype_pool$T_events[e_step$idxs, ]) /
          sum(e_step$w)

        return(c(expected_dist, expected_Tdiff))
      }
      expected_dist = ret[, 1]
      expected_Tdiff = ret[, -1]
      colnames(expected_Tdiff) = NULL
    
    }
    if (any(is.na(expected_Tdiff))) {
      cat("At iteration", iter, "unexpected value for expected time differences\n")
      save(expected_Tdiff,
           file=paste("expected_Tdiff_p", p, "_iter", iter, ".RData", sep=""))
    }
    if (any(is.na(expected_dist))) {
      cat("At iteration", iter, "unexpected value for expected distance\n")
      save(expected_dist,
           file=paste("expected_dist_p", p, "_iter", iter, ".RData", sep=""))
    }
    
    # M step
    eps = mean(expected_dist) / p
    lambdas = 1 / apply(expected_Tdiff, 2, mean)
    if (verbose) cat("M-step - ", iter, "epsilon: ", eps,  "\n")
    if (any(lambdas > max_lambda_val)) {
      idx = which(lambdas > max_lambda_val)
      lambdas[idx] = max_lambda_val
    }
    
    llhood = c(llhood,
               complete_log_likelihood(lambdas, expected_Tdiff, expected_dist,
                                       eps))
    
    if (verbose) {
      cat("Iteration:", iter, "- Log-likelihood: ", llhood[iter], "\n")
    }
    
    if (iter > record_iter) {
      if (verbose) {
        cat("Recording parameters .. \n")
      }
      avg_eps     = avg_eps + eps
      avg_lambdas = avg_lambdas + lambdas
      avg_llhood  = avg_llhood + llhood[iter]
    }
    iter = iter + 1
  }
  avg_eps     = avg_eps / (max_iter - record_iter)
  avg_lambdas = avg_lambdas / (max_iter - record_iter)
  avg_llhood  = avg_llhood / (max_iter - record_iter)
  return(list("lambdas"=lambdas, "eps"=eps, "llhood"=llhood,
              "avg_lambdas"=avg_lambdas, "avg_eps"=avg_eps,
              "avg_llhood"=avg_llhood))
}


complete_log_likelihood_ <- function(lambda, Tdiff, dist, eps, W) {
  .Call('_complete_log_likelihood', PACKAGE = 'mccbn', lambda, Tdiff, dist,
        eps, W)
}

MCEM_hcbn_ <- function(
  lambda, poset, obs, lambda.s, topo.path, L,
  sampling=c('naive', 'add-remove', 'backward'), times=NULL, weights=NULL,
  version=3, perturb.prob=0.3, max.iter=100, burn.in=0.8, max_lambda=1e6,
  thrds=1, verbose=FALSE, seed=NULL) {

  sampling <- match.arg(sampling)
  N <- nrow(obs)
  if (!is.integer(poset)) {
    poset <- matrix(as.integer(poset), nrow=nrow(poset), ncol=ncol(poset))
  }
  if (!is.integer(obs)) {
    obs <- matrix(as.integer(obs), nrow=N, ncol=ncol(obs))
  }
  if (!is.integer(topo.path)) {
    topo.path <- as.integer(topo.path)
  }
  if (is.null(times)) {
    times <- rep(0, N)
    sampling.times.available <- FALSE
  } else {
    sampling.times.available <- TRUE
  }
  if (is.null(weights)) {
    weights <- rep(1, N)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  .Call('_MCEM_hcbn', PACKAGE = 'mccbn', lambda, poset, obs, times,
        lambda.s, topo.path, weights, L, sampling, version, perturb.prob,
        max.iter, burn.in,  max_lambda, sampling.times.available,
        as.integer(thrds), verbose)
}

importance_weight_ <- function(
  genotype, L, poset, lambda, topo.path, eps, time=NULL,
  sampling=c('naive', 'add-remove', 'backward'), version, perturb.prob,
  dist.pool=integer(0), Tdiff.pool=matrix(0), lambda.s=1.0, seed=NULL) {

  sampling = match.arg(sampling)
  if (!is.integer(genotype)) {
    genotype <- as.integer(genotype)
  }
  if (!is.integer(poset)) {
    poset <- matrix(as.integer(poset), nrow=nrow(poset), ncol=ncol(poset))
  }
  if (!is.integer(topo_path)) {
    topo.path <- as.integer(topo_path)
  }
  if (is.null(time)) {
    time <- 0
    sampling.times.available <- FALSE
  } else {
    sampling.times.available <- TRUE
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  .Call('_importance_weight', PACKAGE = 'mccbn', genotype, L, poset,
        lambda, topo.path, eps, time, sampling, version, perturb.prob,
        dist.pool, Tdiff.pool, lambda.s, sampling.times.available)
}
