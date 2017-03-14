hamming_dist <- function(x, y) {
  # For binary genotypes, summing |x-y| 
  dist = sum(abs(x - y))
  #dist = sum(x != y)
  return(dist)
}


possible_genotypes <- function(p) {
  m = 2^p
  genotypes = matrix(0, nrow=m, ncol=p)
  for (i in 1:m) {
    genotypes[i, ] = gen_genotype(p, i)
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


complete_log_likelihood <- function(lambdas, Tdiff, dist, eps) {
  # complete-data log-likelihood or hidden log-likelihood
  p = length(lambdas)
  N = length(dist)
  llhood = N * sum(log(lambdas)) - sum(lambdas * t(Tdiff)) + 
    p * log(eps) * sum(dist) + p * log(1-eps) * sum(p - dist)
  return(llhood)
}


obs_log_likelihood <- function(obs_events, poset, lambdas, lambda_s, eps, 
                               sampling_times=NULL, L=1000, 
                               sampling=c('naive', 'add-remove'), exact=FALSE) {
  
  sampling = match.arg(sampling)
  
  if (exact) {
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
    if (is.null(sampling_times)) {
      prob_Y = apply(obs_events, 1, prob_importance_sampling, L=L, poset=poset,
                     lambdas=lambdas, lambda_s=lambda_s, eps=eps, 
                     sampling_time=NULL, sampling=sampling)
    } else {
      N = nrow(obs_events)
      prob_Y = sapply(1:N, function(i) {
        prob_importance_sampling(genotype=obs_events[i, ], L=L, poset=poset, 
                                 lambdas=lambdas, lambda_s=lambda_s, eps=eps, 
                                 sampling_time=sampling_times[i], sampling=sampling)
      })
    }
    
    llhood = sum(log(prob_Y))
  }
  
  return(llhood)
}


geno_prob_empirical <- function(N, poset, lambdas, lambda_s, genotype, eps) {
  # Compute P(genotype) empirically
  simGenotypes = sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas,
                                  eps=eps)
  return(sum(apply(simGenotypes$obs_events, 1, 
                   function(x, y) all(x == y), y=genotype)) / N)
}


tdiff_empirical <- function(N, poset, lambdas, lambda_s, genotype, eps) {
  # Compute the time differences expirically for a given genotype
  simGenotypes = sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas,
                                  eps=eps)
  idx = which(apply(simGenotypes$obs_events, 1, function(x, y) all(x == y), 
                    y=genotype))
  
  return(apply(simGenotypes$T_events[idx, ], 2, mean))
}


dist_empirical <- function(N, poset, lambdas, lambda_s, genotype, eps) {
  # Compute the average distance between genotype Y (subject to noise) and N
  # possible true genotypes X
  p = ncol(poset)
  simGenotypes = sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas, eps=eps)
  idx = which(apply(simGenotypes$obs_events, 1, function(x, y) all(x == y),  y=genotype))
  
  return(sum(apply(simGenotypes$hidden_genotypes[idx, ], 1, hamming_dist, y=genotype)) / length(idx))
}


get_parents <- function(poset) {
  parents <- apply(poset, 2, function(x) which(x == 1))
  return(parents)
}


get_childreen <- function(poset) {
  childreen <- apply(poset, 1, function(x) which(x == 1))
  return(childreen)
}


add_relations <- function(genotype, childreen) {
  
  # Add 1's:   look for parent nodes which have not occurred, but any of their 
  #            childreen had, and insert a 1 (replace 0 in parental node by 1).
  # genotype:  binary vector indicating whether an event has been observed (1)
  #            or not (0)
  # childreen: computed on the transitively closed poset 
  
  p = length(genotype)
  return(sapply(1:p, 
                function(i, obs, c) 
                  ifelse(obs[i] == 0, ifelse(any(obs[c[[i]]] == 1), 1, obs[i]),
                         obs[i]), obs=genotype, c=childreen))
}


remove_relations <- function(genotype, parents) {
  
  # Remove 1's: look for child nodes which have occurred, but if any of their  
  #             parents have not occurred remove them (replace 1 in child 
  #             node by 0).
  # genotype:   binary vector indicating whether an event has been observed (1)
  #             or not (0)
  # parents:    computed on the transitively closed poset 
  
  p = length(genotype)
  return(sapply(1:p, 
                function(i, obs, p) 
                  ifelse(obs[i] == 1, 
                         ifelse(any(obs[p[[i]]] == 0), 0, obs[i]), obs[i]),
                obs=genotype, p=parents))
}


get_possible_moves <- function(genotype, parents, childreen) {

  # Get the size of the Exit set - mutations that can happen next
  idx = which(genotype == 0)
  if (length(idx) == 0) {
    set_size_add = 0
    idxs_add = NA
  } else {
    idxs_add = idx[sapply(idx, function(i, x, p) all(x[p[[i]]] == 1), x=genotype, 
                          p=parents)]
    set_size_add = length(idxs_add)
  }
  
  # Get the size of the set of mutations that happened last - can be removed
  # and the genotype remains compatible with the poset
  idx = which(genotype == 1)
  if (length(idx) == 0) {
    set_size_remove = 0
    idxs_remove = NA
  } else {
    idxs_remove = idx[sapply(idx, function(i, x, c) all(x[c[[i]]] == 0), 
                             x=genotype, c=childreen)]
    set_size_remove = length(idxs_remove)
  }
  
  return(list("set_size"=c("add"=set_size_add, 
                           "remove"=set_size_remove), 
              "idxs_add"=idxs_add, "idxs_remove"=idxs_remove))
}


perturb <- function(events, parents, childreen, perturb_prob=0.8, 
                    compatible=FALSE) {
  
  # Draw a random number between 0 and 1. If this number is < than perturb_prob,
  # then perturb the genotype by adding or removing an observation, but ensuring
  # observation remains compatible
  
  # ARGUMENTS
  # events:    matrix of genotypes, each row correponds to a vector indicating
  #            whether an event has been observed (1) or not (0)
  # parents:   computed on the transitively closed poset 
  # childreen: computed on the transitively closed poset 
  
  # RETURN VALUES
  # X:                perturbed version of compatible genotypes, remains 
  #                   compatible with the current poset
  # option_set_size:  size of the exit set or set of observations which can be 
  #                   removed. NA is returned when the sample was not perturbed 
  
  L = nrow(events)  # number of samples
  p = ncol(events)  # number of events/mutations
  new_obs = events
  option_set_size = rep(NA, L) 
  
  perturb = ifelse(runif(L, 0, 1) < perturb_prob, 1, 0)
  idxs_perturb = which(perturb == 1)
  
  # When observation Y was already compatible, all rows are equivalent
  if (compatible) {
    # Get indices of mutations that can be either added or removed
    idxs = get_possible_moves(events[1, ], parents, childreen)
  }
  
  for (i in idxs_perturb) {
    
    # Get indices of mutations that can be either added or removed
    if (!compatible) {
      idxs = get_possible_moves(events[i, ], parents, childreen)
    }
    option_set_size[i] = sum(idxs$set_size)
    
    # Choose one index randomly
    idx = sample.int(option_set_size[i], 1)
    mutation_idx = na.omit(c(idxs$idxs_add, idxs$idxs_remove))[idx]
    
    # Check whether move corresponds to an adding move
    add = ifelse(idx <= idxs$set_size["add"], TRUE, FALSE)
    
    new_obs[i, mutation_idx] = ifelse(add, 1, 0)
    
  }
  
  return(list("X" = new_obs, "option_set_size" = option_set_size))
}


draw_samples <- function(genotype, L, parents, childreen, compatible, eps,
                         perturb_prob=0.8) {
  
  # Draw L samples from the proposal
  # two-steps proposal: make compatible and perturb
  # genotype:   binary vector indicating whether an event has been observed (1)
  #             or not (0)
  # L:          number of samples
  # parents:    computed on the transitively closed poset 
  # childreen:  computed on the transitively closed poset 
  # compatible: variable indicating whether the genotype is compatible (1) with 
  #             current poset or not (0)
  # perturb_prob: Genotypes are perturbed in order to learn epsilon. A genotype 
  #             is perturbed, if after drawing a random number, this is <
  #             perturb_prob, otherwise X_start = X
  
  p = length(genotype)  # number of events/mutations
  
  if (compatible) {
    # if genotype is compatible with poset, we can inmediately go to step 2:
    # perturb
    compatible_obs = matrix(genotype, nrow=L, ncol=p, byrow=TRUE)
  } else {
    # if genotype is NOT compatible with poset, we need to generate L compatible
    # versions by adding or removing observations
    compatible_by_adding = add_relations(genotype, childreen)
    compatible_by_removing = remove_relations(genotype, parents)
    
    dist_add = hamming_dist(genotype, compatible_by_adding)
    dist_remove = hamming_dist(genotype, compatible_by_removing)
    add_prob = eps^dist_add / (eps^dist_add + eps^dist_remove)
    add = ifelse(runif(L, 0, 1) < add_prob, 1, 0)
    compatible_obs = matrix(0, L, p)
    
    if (all(add == 1)) {
      compatible_obs = matrix(compatible_by_adding, nrow=L, ncol=p, byrow=TRUE)
    } else if (all(add == 0)) {
      compatible_obs = matrix(compatible_by_removing, nrow=L, ncol=p, byrow=TRUE)
    } else {
      idxs_add = which(add == 1)
      compatible_obs[idxs_add, ] = matrix(compatible_by_adding, nrow=sum(add), 
                                          ncol=p, byrow=TRUE)
      compatible_obs[-idxs_add, ] = matrix(compatible_by_removing, nrow=L-sum(add), 
                                           ncol=p, byrow=TRUE)
    }
  }
  samples = perturb(compatible_obs, parents, childreen, perturb_prob, compatible)
  return(list("samples"=samples, "add_prob"=add_prob))
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


importance_weight <- function(genotype, L, poset, lambdas, lambda_s, eps, 
                              sampling_time=NULL, perturb_prob=0.8,
                              sampling=c('naive', 'add-remove')) {
  
  sampling = match.arg(sampling)
  p = ncol(poset) # number of events/mutations
  
  if (sampling == 'naive') {
    
    # Generate L samples from poset with parameters 'lambdas' and 'lambda_s'. In 
    # particular epsilon is zero (default value) - because the idea is to generate
    # samples of X (underlying true)
    ### **TODO** ### If sampling times available, not considered to generate 
    # samples
    if (!is.null(sampling_time)) {
      warning("Naive proposal doesn't account for sampling times")
    }
    samples = sample_genotypes(L, poset, sampling_param=lambda_s, lambdas=lambdas)
    d = apply(samples$hidden_genotypes, 1, hamming_dist, y=genotype)
    prob_Y_X = eps^d * (1-eps)^(p-d)
    
    return(list("w"=prob_Y_X, "time_differences"=samples$T_events, "dist"=d))
    
  } else if (sampling == 'add-remove') {
    
    poset_trans_closed = trans_closure(poset)
    parents = get_parents(poset_trans_closed)
    childreen = get_childreen(poset_trans_closed)
    compatible = is_compatible(genotype, poset)
    # Generate L samples according to proposal
    return_list = draw_samples(genotype, L, parents, childreen, compatible, eps, 
                               perturb_prob)
    add_prob = return_list$add_prob
    samples = return_list$samples
    
    # Generate mutation times Tj from sample i
    Tdiff = apply(samples$X, 1, sample_mutation_times, poset=poset,
                  lambdas=lambdas, lambda_s=lambda_s, 
                  sampling_time=sampling_time)
    
    log_proposal_X = Tdiff["density", ]
    Tdiff = t(tail(Tdiff, n=-1))
    
    # Hamming distance bewteen samples and genotype
    dist = apply(samples$X, 1, hamming_dist, y=genotype)
    
    # Computing log(Pr(Y|X))
    log_Pr_Y_X = log(eps^dist) + log((1-eps)^(p-dist))
    
    # Computing log(Pr(X))
    log_Pr_X = apply(Tdiff, 1, log_cbn_density, rate=lambdas)
    
    # Computing density of the proposal - for correction
    # proposal : weight accounting for making genotype compatible + weight 
    #            accounting for choosing a possible mutation for perturbing
    #            the sample + weight accounting for sampling times from 
    #            a truncated exponential
    num_options = samples$option_set_size
    stopifnot(all(na.omit(num_options) != 0))
    log_proposal_Y_X = ifelse(is.na(num_options), log(1-perturb_prob), 
                              log(perturb_prob) + log(1/num_options))
    if (!compatible) {
      # TODO: need to keep track which genotype was made compatible by adding or 
      # removing -- e.g. instead of returnig add_prob, return compatible_prob
      log_proposal_Y_X = log_proposal_Y_X + log(add_prob) 
    } 
    
    log_proposal = log_proposal_Y_X + log_proposal_X
    
    importance_weight = exp(log_Pr_X + log_Pr_Y_X - log_proposal)
    return(list("w"=importance_weight, "time_differences"=Tdiff, "dist"=dist))
  }

}


prob_importance_sampling <- function(genotype, L, poset, lambdas, lambda_s, 
                                     eps, sampling_time=NULL, 
                                     sampling=c('naive', 'add-remove')) {
  
  # Compute Pr(genotype) using Monte Carlo sampling
  # L        number of samples
  # sampling type of sampling. OPTIONS: 'naive' - generating occurrence times  
  #          according to lambdas and from them genotypes; 'add-remove' - 
  #          generating genotypes from observed genotypes using a two-steps 
  #          proposal. First, making genotypes compatible with the poset by
  #          either adding or removing observations. Second, perturbing this 
  #          version by adding or removing observations. 
  # NOTE: If 'naive' sampling is employed, sampling times are not used. 
  
  sampling = match.arg(sampling)
  probs = importance_weight(genotype, L, poset, lambdas, lambda_s, eps, 
                            sampling_time=sampling_time, sampling=sampling)
  probs = probs$w

  return(sum(probs) / L)
}


tdiff_importance_sampling <- function(genotype, L, poset, lambdas, lambda_s, 
                                      eps, sampling_time=NULL, 
                                      sampling=c('naive', 'add-remove')) {
  
  # Compute expected time differences for a given genotype using Monte Carlo 
  # sampling
  sampling = match.arg(sampling)
  importance_weights = importance_weight(genotype, L, poset, lambdas, lambda_s,
                                         eps, sampling_time=sampling_time, 
                                         sampling=sampling)
  return(colSums(importance_weights$w * importance_weights$time_differences) / 
           sum(importance_weights$w))

}


dist_importance_sampling <- function(genotype, L, poset, lambdas, lambda_s, 
                                     eps, sampling_time=NULL, 
                                     sampling=c('naive', 'add-remove')) {
  
  # Compute expected hamming distance between a given observations (genotype) 
  # and the underlying/true genotype using Monte Carlo sampling
  sampling = match.arg(sampling)
  importance_weights = importance_weight(genotype, L, poset, lambdas, lambda_s,
                                         eps, sampling_time=sampling_time, 
                                         sampling=sampling)
  return(sum(importance_weights$dist * importance_weights$w) / 
           sum(importance_weights$w))
}


prob_empirical_vs_sampling <- function(events, L, rep=NULL, one_genotype=FALSE, 
                                  sampling_times=NULL, 
                                  sampling=c('naive', 'add-remove'), outdir,
                                  outname="", binwidth=0.01) {
  
  # events         matrix of observations or single genotype, if 'one_genotype' 
  #                is TRUE 
  # L              number of samples
  # rep            if only one genotype is provided, number of repetitions
  # one_genotype   boolean variable indicating if only one genotype is provided
  # sampling_times vector of sampling times per event. 
  # sampling       type of sampling. OPTIONS: 'naive' or 'add-remove'.
  # NOTE: If 'naive' sampling is employed, sampling times are not used. 
  
  sampling = match.arg(sampling)
  
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
  
  prob_empirical = numeric(N)
  prob_sampling  = numeric(N)
  
  for (i in 1:N) {
    if (!one_genotype) {
      genotype = simulated_obs$obs_events[i, ]
      sampling_time = sampling_times[i]
    }
    prob_empirical[i] = geno_prob_empirical(N=100000, poset, lambdas, lambda_s, 
                                            genotype, eps=eps)
    prob_sampling[i] = prob_importance_sampling(genotype, L=L, poset, lambdas, 
                                                lambda_s, eps=eps, 
                                                sampling_time=sampling_time,
                                                sampling=sampling)
  }
  
  outdir = file.path(outdir, paste("L", L, "_", sampling, sep=""))
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  outname = file.path(outdir, paste("probability_Y_empirical_vs_sampling", 
                                    outname, ".pdf", sep=""))
  
  if (one_genotype) { 
    xlab = expression(P(Y))
    ylab = ""
  } else {
    xlab = expression(P[empirical](Y))
    ylab = expression(widehat(P)(Y))
  }
  
  pl_empirical_vs_sampling(empirical=prob_empirical, sampling=prob_sampling, 
                           xlab=xlab, ylab=ylab, N=N, one_genotype=one_genotype,
                           outname=outname, binwidth=binwidth)
  
  return(list("empirical" = prob_empirical, "sampling" = prob_sampling))
}


tdiff_empirical_vs_sampling <- function(events, L, rep=NULL, one_genotype=FALSE, 
                                       sampling_times=NULL, 
                                       sampling=c('naive', 'add-remove'), outdir,
                                       outname="", binwidth=0.01) {
  
  # events         matrix of observations or single genotype, if 'one_genotype' 
  #                is TRUE 
  # L              number of samples
  # rep            if only one genotype is provided, number of repetitions
  # one_genotype   boolean variable indicating if only one genotype is provided
  # sampling_times vector of sampling times per event. 
  # sampling       type of sampling. OPTIONS: 'naive' or 'add-remove'.
  # NOTE: If 'naive' sampling is employed, sampling times are not used. 
  
  sampling = match.arg(sampling)
  
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
    time_diff_empirical[i, ] = tdiff_empirical(N=100000, poset, lambdas, lambda_s,
                                               genotype, eps=eps)
    time_diff_sampling[i, ] = tdiff_importance_sampling(genotype, L=L, poset, 
                                                        lambdas, lambda_s, eps,
                                                        sampling_time=sampling_time, 
                                                        sampling=sampling)
  }
  
  outdir = file.path(outdir, paste("L", L, "_", sampling, sep=""))
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  
  for (j in 1:p) {
    pl_name = file.path(outdir, paste("time_diff_empirical_vs_sampling", 
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
  
  return(list("empirical" = time_diff_empirical, "sampling" = time_diff_sampling))
}


dist_empirical_vs_sampling <- function(events, L, rep=NULL, one_genotype=FALSE, 
                                       sampling_times=NULL, 
                                       sampling=c('naive', 'add-remove'), outdir,
                                       outname="", binwidth=0.01) {
  
  # events         matrix of observations or single genotype, if 'one_genotype' 
  #                is TRUE 
  # L              number of samples
  # rep            if only one genotype is provided, number of repetitions
  # one_genotype   boolean variable indicating if only one genotype is provided
  # sampling_times vector of sampling times per event. 
  # sampling       type of sampling. OPTIONS: 'naive' or 'add-remove'.
  # NOTE: If 'naive' sampling is employed, sampling times are not used. 
  
  sampling = match.arg(sampling)
  
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
    d_empirical[i] = dist_empirical(N=100000, poset, lambdas, lambda_s, genotype,
                                    eps=eps)
    d_sampling[i] = dist_importance_sampling(genotype, L=L, poset, lambdas, 
                                             lambda_s, eps, 
                                             sampling_time=sampling_time,
                                             sampling=sampling)
  }
  
  outdir = file.path(outdir, paste("L", L, "_", sampling, sep=""))
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  outname = file.path(outdir, paste("hamming_dist_empirical_vs_sampling", 
                                    outname, ".pdf", sep=""))
  
  if (one_genotype) { 
    xlab = expression(d(X,Y))
    ylab = ""
  } else {
    xlab = expression(d[empirical](X,Y))
    ylab = expression(widehat(d)(X,Y))
  }
  
  pl_empirical_vs_sampling(empirical=d_empirical, sampling=d_sampling, xlab=xlab, 
                           ylab=ylab, N=N, one_genotype=one_genotype, 
                           outname=outname, binwidth=binwidth)
  
  return(list("empirical" = d_empirical, "sampling" = d_sampling))
}


pl_empirical_vs_sampling <- function(empirical, sampling, xlab, ylab, 
                                     N=NULL, one_genotype=FALSE, outname="", 
                                     binwidth=0.01) {
  
  if (one_genotype) {
    df = data.frame(x=c(empirical, sampling), 
                    method=c(rep("empirical", N), rep("sampling", N)))
    pl = ggplot(df, aes(x = x, fill = method))
    pl = pl + geom_histogram(binwidth=binwidth, alpha=0.5, position="identity") + 
      labs(x=xlab, y="Frequency") + 
      theme_bw() + theme(text=element_text(size=14))
    
    ggsave(outname, pl, width=4, height=2.5)
    
  } else {
    df = data.frame(x = empirical, y = sampling)
    pl = ggplot(df, aes(x = x, y = y))
    pl = pl + geom_point() + 
      geom_abline(intercept = 0, slope = 1, colour="blue") + 
      labs(x = xlab, y = ylab) + 
      theme_bw() + theme(text=element_text(size=14))
    
    ggsave(outname, pl, width=3, height=2)
  }
  
}


MCEM_hcbn <- function(poset, obs_events, sampling_times=NULL, lambda_s=1.0,   
                      max_iter=100, burn_in=0.8, L=100, 
                      sampling=c('naive', 'add-remove'), max_lambda_val=1e6, 
                      perturb_prob=0.8, parallel=TRUE, verbose=TRUE)  {
  
  # obs_events     matrix of observations, each row correponds to a vector 
  #                indicating whether an event has been observed (1) or not (0)
  # sampling_times vector of sampling times per observation/genotype 
  # L              number of samples to be drawn from the proposal
  # sampling       type of sampling. OPTIONS: 'naive' or 'add-remove'.
  # perturb_prob   Genotypes are perturbed in order to learn epsilon. A genotype 
  #                is perturbed, if after drawing a random number, this is <
  #                perturb_prob, otherwise X_start = X. Option is used if  
  #                sampling is set to 'add-remove'.
  # parallel       boolean variable indicating whether sampling should be
  #                executed sequentially (0) or in parallel (1). Option is used if
  #                sampling is set to 'naive'.
  # NOTE: If 'naive' sampling is employed, sampling times are not used.
  
  sampling = match.arg(sampling)
  
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
  
  # initialize parameters
  eps = runif(1, 0, 0.5)
  lambdas = initialize_lambda(obs_events=obs_events, 
                              average_sampling_times=avg_sampling_t, 
                              poset=poset, verbose=verbose)
  
  # NOTE: function initialize_lambda should be modify to account for cases where
  #       sum(obs_events[parents_i == 1, i]) is zero. 
  if (any(lambdas == 0))  {
    idxs = which(lambdas == 0)
    lambdas[idxs] = 1e-6
  }
  
  avg_eps = 0 
  avg_lambdas = numeric(p)
  avg_llhood = 0
  
  record_iter = max(as.integer(burn_in * max_iter), 1)
  
  if (sampling == 'add-remove') {
    poset_trans_closed = trans_closure(poset)
    parents = get_parents(poset_trans_closed)
    childreen = get_childreen(poset_trans_closed)
    # TODO: at the moment, compatibility test is performed for each genotype
    #idx_compatible = apply(obs_events, 1, is_compatible, poset=poset)
  }
  
  for(iter in 1:max_iter) {
    
    # E step
    # Compute for each event j, expected_Tdiff[j] = E[T_j - max_{u \in pa(j)} T_u | Y]
    # Compute expected_dist
    
    if (sampling == 'naive') {
      
      if (!parallel) {
        expected_dist = numeric(N)
        expected_Tdiff = matrix(0, nrow=N, ncol=p)
        for (i in 1:N) {
          # Draw L samples from poset with parameters 'lambdas' and 'lambda_s'
          e_step = importance_weight(genotype=obs_events[i, ], L=L, poset=poset,
                                     lambdas=lambdas, lambda_s=lambda_s, eps=eps,
                                     sampling_time=NULL, sampling='naive')
          # Conditional expectation of the sufficient statistic d(X, Y)
          expected_dist[i] = sum(e_step$w * e_step$dist) / sum(e_step$w)
          
          # Contitional expectation of the sufficient statistic Z_j
          # Z_j = t_j - max_{u \in pa(j)} t_u
          expected_Tdiff[i, ] = colSums(e_step$w * e_step$time_differences) /
            sum(e_step$w)
        }
      } else {
        ret = foreach(i=1:N, .combine='rbind') %dopar% {
          
          # Draw L samples from poset with parameters 'lambdas' and 'lambda_s'
          e_step = importance_weight(genotype=obs_events[i, ], L=L, poset=poset,
                                     lambdas=lambdas, lambda_s=lambda_s, eps=eps,
                                     sampling_time=NULL, sampling='naive')
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
      
      ret = foreach(i=1:N, .combine='rbind') %dopar% {
        
        # In each iteration and for each observation, draw L samples from proposal
        # two-steps proposal: make compatible and perturb
        # 1. Make genotypes compatible by adding or removing 1's
        # 2. Perturbe version of previous genotypes, but ensure it remains
        #    compatible with current poset
        e_step = importance_weight(genotype=obs_events[i, ], L=L, poset=poset,
                                   lambdas=lambdas, lambda_s=lambda_s, eps=eps,
                                   sampling_time=sampling_times[i], 
                                   perturb_prob=perturb_prob, sampling='add-remove')
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
    
    if (any(is.na(expected_Tdiff))) {
      stop("Unexpected value for expected time differences")
    }
    if (any(is.na(expected_dist))) {
      stop("Unexpected value for expected time differences")
    }
    
    # M step
    eps = mean(expected_dist) / p
    lambdas = 1 / apply(expected_Tdiff, 2, mean)
    
    if (any(lambdas > max_lambda_val)) {
      idx = which(lambdas > max_lambda_val)
      lambdas[idx] = max_lambda_val
    }
    
    llhood = complete_log_likelihood(lambdas, expected_Tdiff, expected_dist, eps)
    
    if (verbose) {
      cat("Iteration:", iter, "- Log-likelihood: ", llhood, "\n")
    }
    
    if (iter > record_iter) {
      if (verbose) {
        cat("Recording parameters .. \n")
      }
      avg_eps     = avg_eps + eps
      avg_lambdas = avg_lambdas + lambdas
      avg_llhood  = avg_llhood + llhood
    }
  }
  avg_eps     = avg_eps / (max_iter - record_iter)
  avg_lambdas = avg_lambdas / (max_iter - record_iter)
  avg_llhood  = avg_llhood / (max_iter - record_iter)
  return(list("lambdas"=lambdas, "eps"=eps, "llhood"=llhood, "avg_lambdas"=avg_lambdas, 
              "avg_eps"=avg_eps, "avg_llhood"=avg_llhood))
}


