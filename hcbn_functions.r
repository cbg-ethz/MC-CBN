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
  p = length(lambdas)
  n = length(dist)
  llhood = N * sum(log(lambdas)) - sum(lambdas * t(Tdiff)) + 
    p * log(eps) * sum(dist) + p * log(1-eps) * sum(p - dist)
  return(llhood)
}


obs_log_likelihood <- function(obs_events, poset, lambdas, lambda_s, eps, 
                               L=1000, exact=FALSE) {
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
    prob_Y = apply(obs_events, 1, prob_imp, L=L, poset=poset, lambdas=lambdas, 
                   lambda_s=lambda_s, eps=eps)
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





tdiff_imp <- function(L, poset, lambdas, lambda_s, genotype, eps) {
  # Compute expected time differences for a given genotype using Monte Carlo sampling
  p = ncol(poset)
  # Generate L samples from poset with parameters 'lambdas' and 'lambda_s'. In 
  # particular epsilon is zero (default value)
  simGenotypes = sample_genotypes(L, poset, sampling_param=lambda_s, lambdas=lambdas)
  dist = apply(simGenotypes$hidden_genotypes, 1, hamming_dist, y=genotype)
  probs = eps^dist * (1-eps)^(p-dist)
  
  return(apply(simGenotypes$T_events, 2, function(x) sum(x * probs)) / sum(probs))
}


dist_imp <- function(L, poset, lambdas, lambda_s, genotype, eps) {
  # Compute expected hamming distance between a given observations (genotype) and the
  # underlying/true genotype using Monte Carlo sampling
  p = ncol(poset)
  simGenotypes = sample_genotypes(L, poset, sampling_param=lambda_s, lambdas=lambdas)
  dist = apply(simGenotypes$hidden_genotypes, 1, hamming_dist, y=genotype)
  probs = eps^dist * (1-eps)^(p-dist)
  return(sum(dist * probs) / sum(probs))
}


MCMC_hcbn <- function(poset, obs_events, sampling_times=NULL, lambda_s=1.0, max_iter=100,  
                      burn_in=0.2, L=100, max_lambda_val=1e6, verbose=TRUE)  {
  
  p = ncol(poset)       # number of events/mutations
  N = nrow(obs_events)  # number of observations/genotypes
  
  if (is.null(sampling_times)) {
    sampling_times = rep(0, nrow(obs_events))
    avg_sampling_t = lambda_s
  } else {
    avg_sampling_t = mean(sampling_times)
  }
  
  # initialize parameters
  eps = runif(1, 0, 0.5)
  lambdas = initialize_lambda(obs_events=obs_events, average_sampling_times=avg_sampling_t,
                              poset=poset, verbose=verbose)
  
  avg_eps = 0 
  avg_lambdas = numeric(p)
  avg_llhood = 0
  
  record_iter = max(as.integer((1-burn_in) * max_iter), 1)
  
  for(iter in 1:max_iter) {
    
    # E step
    # Compute for each event j, expected_Tdiff[j] = E[T_j - max_{u \in pa(j)} T_u | Y]
    # Compute expected_dist
    expected_dist = numeric(N)
    expected_Tdiff = matrix(0, nrow=N, ncol=p)
    for(i in 1:N) {
      # Draw L samples from poset with parameters "lambdas" 
      simGenotypes = sample_genotypes(L, poset=poset, sampling_param=lambda_s,
                                      lambdas=lambdas)
      
      # Conditional expectation of the sufficient statistic d(X, Y)
      dist = apply(simGenotypes$obs_events, 1, hamming_dist, y=obs_events[i, ])
      pr_Y_X = eps^dist * (1-eps)^(p-dist) 
      expected_dist[i] = sum(pr_Y_X * dist) / sum(pr_Y_X)
      
      # Contitional expectation of the sufficient statistic Z_j
      # Z_j = t_j - max_{u \in pa(j)} t_u
      expected_Tdiff[i, ] = apply(simGenotypes$T_events, 2,
                                  function(x) sum(x * pr_Y_X)) / sum(pr_Y_X)
      
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


get_next_mutation <- function(genotype, parents) {
  # Look for mutations that haven't occurred but its parents have. These 
  # mutations can occur next.
  idx = which(genotype == 0)
  if (length(idx) == 0) {
    exit_size = 0
    mutation_idx = NA
  } else {
    exit = idx[sapply(idx, function(i, x, p) all(x[p[[i]]] == 1), x=genotype, 
                      p=parents)]
    exit_size = length(exit)
    if (exit_size > 1) {
      mutation_idx = sample(exit, 1)
    } else {
      mutation_idx = exit
    }
  }
  return(c("exit_size"=exit_size, "mutation_idx"=mutation_idx))
}


perturb_add <- function(events, idxs) {
  # add 1: add an observation from the Exit(obs) - mutations that can occur next
  L = nrow(events)  # number of samples 
  if (is.null(L)) {
    # only one observation is passed as 'events'
    new_obs = events
    new_obs[idxs] = 1
  } else {
    new_obs = t(sapply(1:L, function(i, x, idxs) { 
      x[i, idxs[i]] = 1
      return(x[i, ])
    }, x=events, idxs=idxs)) 
  }
  return(new_obs)
}


get_previous_mutation <- function(genotype, childreen) {
  # Look for mutations that have occurred but any of its childreen haven't. These
  # mutation are candidates to be removed
  idx = which(genotype == 1)
  if (length(idx) == 0) {
    set_size = 0
    mutation_idx = NA
  } else {
    mutation_idx = idx[sapply(idx, function(i, x, c) all(x[c[[i]]] == 0), 
                              x=genotype, c=childreen)]
    set_size = length(mutation_idx)
    if (set_size > 1) {
      mutation_idx = sample(mutation_idx, 1)
    } 
  }
  return(c("set_size"=set_size, "mutation_idx"=mutation_idx))
}


perturb_remove <- function(events, idxs) {
  # remove 1: remove an observation whose childreen have not occurred
  L = nrow(events)  # number of samples
  if (is.null(L)) {
    new_obs = events
    new_obs[idxs] = 0
  } else {
    new_obs = t(sapply(1:L, function(i, x, idxs) { 
      x[i, idxs[i]] = 0
      return(x[i, ])
    }, x=events, idxs=idxs))
  }
  return(new_obs)
}


perturb <- function(events, parents, childreen, perturb_prob=0.8, add_prob=0.5) {
  
  # Draw a random number between 0 and 1. If this number is < than perturb_prob,
  # then perturb the genotype by adding or removing an observation, but ensuring
  # observation remains compatible
  # ARGUMENTS
  # events:    matrix of genotypes, each row correponds to a vector indicating
  #            whether an event has been observed (1) or not (0)
  # parents:   computed on the transitively closed poset 
  # childreen: computed on the transitively closed poset 
  # RETURN VALUES
  # X:         perturbed version of compatible genotypes, remains compatible
  # set_size:  size of the exit set or set of observations which can be removed
  #            NA is returned when the sample was not perturbed. 
  
  L = nrow(events)  # number of samples
  p = ncol(events)  # number of events/mutations
  new_obs = matrix(0, L, p)
  set_size = rep(NA, L) #NULL #rep(NULL, L)
  
  perturb = ifelse(runif(L, 0, 1) < perturb_prob, 1, 0)
  idxs_perturb = which(perturb == 1)

  if (length(idxs_perturb) == 0) {
    new_obs = events
  } else {
    new_obs[-idxs_perturb, ] = events[-idxs_perturb, ]
  
    # Draw a random number between 0 and 1. If this number is > than add_prob, then 
    # add an observation
    add = ifelse(runif(sum(perturb), 0, 1) > add_prob, 1, 0)
    idxs_add = which(add == 1)
  
    # If length(idxs_perturb) is one, apply returns an error
    events_mat = matrix(events[idxs_perturb, ], nrow=length(idxs_perturb), ncol=p)
  
    if (all(add == 1)) {
    
      # If all samples are perturbed by adding mutations
      add_mutations = apply(events_mat, 1, get_next_mutation, parents=parents)
      new_obs[idxs_perturb, ] = perturb_add(events_mat, 
                                            add_mutations["mutation_idx", ])
      set_size[idxs_perturb] = add_mutations["exit_size", ]
    
    } else if (all(add == 0)) {
    
      # If all samples are perturbed by removing mutations
      remove_mutations = apply(events_mat, 1, get_previous_mutation, 
                               childreen=childreen)
      new_obs[idxs_perturb, ] = perturb_remove(events_mat, 
                                               remove_mutations["mutation_idx", ])
      set_size[idxs_perturb] = remove_mutations["set_size", ]
      
    } else {

      if (length(idxs_add) > 1) {
        # If multiple samples are perturbed by adding mutations
        add_mutations = apply(events_mat[idxs_add, ], 1, get_next_mutation, 
                              parents=parents)
        new_obs[idxs_perturb[idxs_add], ] = perturb_add(events_mat[idxs_add, ], 
                                                        add_mutations["mutation_idx", ])
        set_size[idxs_perturb[idxs_add]] = add_mutations["exit_size", ]
        
      } else {
        # If only one sample is perturbed by adding mutations
        add_mutations = get_next_mutation(events_mat[idxs_add, ], parents=parents)
        new_obs[idxs_perturb[idxs_add], ] = perturb_add(events_mat[idxs_add, ], 
                                                        add_mutations["mutation_idx"])
        set_size[idxs_perturb[idxs_add]] = add_mutations["exit_size"]
        
      }
    
      if ((sum(perturb) - length(idxs_add)) > 1) {
        # If multiple samples are perturbed by removing mutations
        remove_mutations = apply(events_mat[-idxs_add, ], 1, get_previous_mutation, 
                                 childreen=childreen)
        new_obs[idxs_perturb[-idxs_add], ] = perturb_remove(events_mat[-idxs_add, ], 
                                                            remove_mutations["mutation_idx", ])
        set_size[idxs_perturb[-idxs_add]] = remove_mutations["set_size", ]
        
      } else {
        # If only one sample is perturbed by removing mutations
        remove_mutations = get_previous_mutation(events_mat[-idxs_add, ], 
                                                 childreen=childreen)
        new_obs[idxs_perturb[-idxs_add], ] = perturb_remove(events_mat[-idxs_add, ], 
                                                            remove_mutations["mutation_idx"])
        set_size[idxs_perturb[-idxs_add]] = remove_mutations["set_size"]
 
      }
    }
  }
  return(list("X" = new_obs, "set_size" = set_size))
}


draw_samples <- function(genotype, L, parents, childreen, compatible, 
                            add_prob1=0.5, perturb_prob=0.8, add_prob2=0.5) {
  
  # Draw L samples from the proposal
  # two-steps proposal: make compatible and perturb
  # genotype:   binary vector indicating whether an event has been observed (1)
  #             or not (0)
  # L:          number of samples
  # parents:    computed on the transitively closed poset 
  # childreen:  computed on the transitively closed poset 
  # compatible: variable indicating whether the genotype is compatible (1) with 
  #             current poset or not (0)
  # add_prob1:  Incompatible genotypes are made compatible by adding or removing
  #             observations. One of the two operations is chosen by drawing
  #             a random number, and if this number is > add_prob1, then 
  #             observations are added.
  # add_prob2:  Genotypes are pertubed by adding or removing an observation. One 
  #             of the two operations is chosen by drawing a random number, and 
  #             if this number is > add_prob2, then an observation is added. 
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
    add = ifelse(runif(L, 0, 1) > add_prob1, 1, 0)
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
  samples = perturb(compatible_obs, parents, childreen, perturb_prob, add_prob2)
  return(samples)
}


rtexp <- function(x, rate) {
  rand_num = runif(1, 0, 1)
  return(-log(1 - rand_num * (1 - exp(-rate * x))) / rate)
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
  dens = sum(log(lambdas)) - sum(lambdas * Tdiff)
  return(dens)
}


importance_weight <- function(genotype, L, poset, lambdas, lambda_s, eps, 
                              sampling_time, add_prob1=0.5, perturb_prob=0.8, 
                              add_prob2=0.5) {
  p = ncol(poset) # number of events/mutations
  
  poset_trans_closed = trans_closure(poset)
  parents = get_parents(poset_trans_closed)
  childreen = get_childreen(poset_trans_closed)
  compatible = is_compatible(genotype, poset)
  # Generate L samples from poset with parameters 'lambdas' and 'lambda_s' 
  # according to proposal
  samples = draw_samples(genotype, L, parents, childreen, compatible)
  
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
  #            accounting for perturbing by adding or removing + weight  
  #            accounting for choosing a possible event + weight accounting
  #            for sampling times
  set_size = samples$set_size
  set_size[which(set_size == 0)] = 1
  if (compatible) {
    log_proposal_Y_X = ifelse(is.na(set_size), log(1-perturb_prob), 
                              log(perturb_prob) + log(add_prob2) + 
                                log(1/set_size))
    
  } else {
    log_proposal_Y_X = log(add_prob1) + 
      ifelse(is.na(set_size), log(1-perturb_prob), 
             log(perturb_prob) + log(add_prob2) +  log(1/set_size))
  }
  log_proposal = log_proposal_Y_X + log_proposal_X
  
  importance_weight = exp(log_Pr_X + log_Pr_Y_X - log_proposal)
  return(list("w"=importance_weight, "time_differences"=Tdiff))
}


## **TODO** ##
# merge prob_imp and prob_importance_sampling. E.g., include argument 'sampling'
# with options 'random' and 'ar' (stands for the type of moves, add and remove)
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
  
  sampling = match.arg(sampling)
  
  if (sampling == 'naive') {
    
    # Generate L samples from poset with parameters 'lambdas' and 'lambda_s'. In 
    # particular epsilon is zero (default value) - because the idea is to generate
    # samples of X (underlying true)
    samples = sample_genotypes(L, poset, sampling_param=lambda_s, lambdas=lambdas)
    p = ncol(poset)
    d = apply(samples$hidden_genotypes, 1, hamming_dist, y=genotype)
    probs = eps^d * (1-eps)^(p-d)
    
  } else if (sampling == 'add-remove') {
    
    probs = importance_weight(genotype, L, poset, lambdas, lambda_s, eps, 
                              sampling_time)
    probs = probs$w
    
  }
  return(sum(probs) / L)
}


tdiff_importance_sampling <- function(genotype, L, poset, lambdas, lambda_s, eps,
                                      sampling_time=NULL) {
  
  # Compute expected time differences for a given genotype using Monte Carlo 
  # sampling
  
  importance_weights = importance_weight(genotype, L, poset, lambdas, lambda_s, eps, 
                                         sampling_time)
  
  return(colSums(importance_weights$w * importance_weights$time_differences) / 
           sum(importance_weights$w))
}