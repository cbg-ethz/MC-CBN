compute_geno_indexes <- function(obs_events, poset) {
  if(is.matrix(obs_events) == FALSE) {
    obs_events = as.matrix(obs_events)
  }
    
    
  G = compatible_genotypes_with_matrix(poset)
  
  geno_indexes = rep(0, nrow(obs_events))
  for(i in 1:nrow(obs_events)) {
    genotype = obs_events[i, ]
    geno_index = which(apply(G, 1, function(x) { all(genotype==x)} ))
    geno_indexes[i] = ifelse(length(geno_index) == 0, NA, geno_index)
  }
  geno_indexes
}


log_observed_geno_probs <- function(lambdas, obs_events, dtimes, G, geno_indexes) {
  times = dtimes$times
  t_indexes = dtimes$t_indexes
  all_times = dtimes$all_times
  if(length(lambdas) == 1) {
    log_observed_geno_probs = log( abs(obs_events - exp(-lambdas*all_times) ) )
    log_observed_geno_probs[is.infinite(log_observed_geno_probs)] = -50
    return(  log_observed_geno_probs )
  }
  
  if(any(is.na(geno_indexes) ) ) {
    
    print("Some genotypes are not compatible with the poset. Hence, the log likelihood is a very small number (minus infinity)")
    # -Inf does not work with numerical optimization. Hence, -50 multiply by number of genotypes is used for very small likelihood 
    return(rep(-50, length(all_times)))
    
  }
  
  # computing the genotype probabilities for all time points (and all possible genotypes). The return value
  # is a matrix 
  geno_probs = probs(lambdas, G, times)
  
  # computing the genotype probabilities for all time points only for  the observed genotypes. The return value
  # is a vector
  log_observed_geno_probs = log( sapply(1:nrow(obs_events), function(i) { geno_probs[t_indexes[i], geno_indexes[i]] }) )
  
#   print(log_observed_geno_probs)
  # Handling a special case: if the probability for a compatible genotype is zero (for any reason)
  log_observed_geno_probs[is.infinite(log_observed_geno_probs)] = -50
  log_observed_geno_probs   
}


likelihood <- function(lambdas, obs_events, dtimes, G, geno_indexes, ll_control=list(negate=FALSE)) {
  negate = ll_control$negate
  C = ifelse(negate, -1, 1)
  
  C*sum(log_observed_geno_probs(lambdas, obs_events, dtimes, G, geno_indexes))  
}


#####################################  add error modeling
log_observed_geno_probs_with_eps<- function(lambda, eps, obs_events, dtimes, G) {
  obs_events = as.matrix(obs_events, ncol=length(lambda))
  times = dtimes$times
  t_indexes = dtimes$t_indexes
  all_times = dtimes$all_times
  
  # computing the genotype probabilities for all time points (and all possible genotypes). The return value
  # is a matrix 
  geno_probs = probs(lambda, G, times)
  #   range(geno_probs)
  prob_Y = rep(0, nrow(obs_events))
  for( i in 1:nrow(obs_events)) {
    Y = obs_events[i, ]
    
    P_Y_X = apply(G, 1, Y=Y, eps=eps, function(x, Y, eps) { 
      prob_hamming_distance(x, Y, eps)
    })
    
    prob_Y[i] = sum(P_Y_X * geno_probs[t_indexes[i], ] )
  }
  
  
  log_observed_geno_probs = log(prob_Y)
  
  # Handling a special case: if the probability for a genotype is zero (for any reason)
  log_observed_geno_probs[is.infinite(log_observed_geno_probs)] = -50
  
  log_observed_geno_probs
}

likelihood_with_eps <- function(lambda, obs_events, dtimes, G, geno_indexes, ll_control=list(negate=FALSE, eps=0.25)) {
  
  negate = ll_control$negate
  eps = ll_control$eps
  
  log_observed_geno_probs = log_observed_geno_probs_with_eps(lambda, eps, obs_events, dtimes, G)
  C = ifelse(negate, -1, 1)
  C * sum(log_observed_geno_probs)
}



# Y: observed genotype
# X: true genotype
# eps: per-locus error rate
prob_hamming_distance <- function(X, Y, eps) {
  p = length(X)
  d = sum( X != Y )
  (eps^d) * (1-eps)^(p-d)
}


# 
# # C : nr of compatible genotypes
# loglikelihood_cbn <- function(poset, compatible_ll, C, N){
#   p = ncol(poset)
#   lattice_size = orderIdeals(poset)$lattice_size
#   
#   # fraction of the genotypes that are compatible with the given poset
#   alpha = C /N
#   
#   q_e = 1/(2^p - lattice_size  )
#   
#   incompatible_ll = 0.0
#   if((N-C) > 0)
#     incompatible_ll = (N-C)* (log(1-alpha) + log(q_e))
#   
#   compatible_ll = compatible_ll +  C * log(alpha)
#   
#   list(ll=incompatible_ll + compatible_ll, incompatible_ll=incompatible_ll,  compatible_ll=compatible_ll, alpha=alpha)  
# }


compute_q_e <- function(poset) {
  p=ncol(poset)
  if(ncol(poset) < 40 || sum(poset) < 8) {
    log_lattice_size = genoLatticeSize_fast(poset)
    log_q_e = -logsumexp( c(p * log(2), log_lattice_size), c(1, -1) )$logR
  } else {
    log_q_e = -  ncol(poset) * log(2)
  }
  #max(log_q_e, -  ncol(poset)/3 * log(2) )
  log_q_e
}

incompatible_loglike <- function(poset, obs_events, times, weights, compatible_geno) {
  p = ncol(poset)

#   log_lattice_size = genoLatticeSize_fast(poset)
  compatible_geno = compatible_genotypes(obs_events, poset)
  
  # C: number of compatible genotypes
  C = sum(weights[compatible_geno$compatible_indexes])
  N = sum(weights)
  # I: number of incompatible genotypes
  I = N - C
  
  # fraction of the data that are compatible with the given poset
  alpha = C / N
  
  log_q_e = compute_q_e(poset)
  
  incompatible_ll = 0.0
  if( I > 0 ){
    incompatible_ll = I * (log(1-alpha) + log_q_e)
  }
  
  list(ll=incompatible_ll,  alpha=alpha) 
}


#' loglike_mixture_model
#' @export
loglike_mixture_model <- function(poset, lambda, obs_events, times, weights, control = list(ll_method="importance"), compatible_geno, incomp_loglike){
#   p = ncol(poset)
#   
# #   lattice_size = orderIdeals(poset)$lattice_size
#   log_lattice_size = genoLatticeSize(poset)
# #   compatible_geno = compatible_genotypes(obs_events, poset)
# 
#     # C: number of compatible genotypes
#     C = sum(weights[compatible_geno$compatible_indexes])
#     N = sum(weights)
#     # I: number of incompatible genotypes
#     I = N - C
#     
#     # fraction of the data that are compatible with the given poset
#     alpha = C / N
#       
# 
#   log_q_e = -logsumexp( c(p * log(2), log_lattice_size), c(1, -1) )$logR
#   
#   
#   incompatible_ll = 0.0
#   if( I > 0 ){
#     incompatible_ll = I * (log(1-alpha) + log_q_e)
#   }
#     
      
  incompatible_ll = incomp_loglike$ll
  C = sum(weights[compatible_geno$compatible_indexes])
  alpha = incomp_loglike$alpha
  
  compatible_ll = 0.0
  if( C > 0 )
  {
    genotypes = obs_events[compatible_geno$compatible_indexes, , drop=FALSE ]
    
    times = times[compatible_geno$compatible_indexes]
    weights = weights[compatible_geno$compatible_indexes]
    
    if(control$ll_method == "importance") {
       tmp = loglike_importance_sampling(poset, lambda, genotypes, times, control$nrOfSamples, weights, with_eps=FALSE, eps=NA)
       compatible_ll = tmp$approx_loglike
    } else {
      dtimes = discretize_times(times, control$D)
      compatible_ll = .likelihood_decomp(poset, lambda, genotypes, dtimes)
    }
    
    compatible_ll = compatible_ll + C * log(alpha)
  }
  
#   
#   ################### temporary to remove
#   print("todo")
#   tmp2 = incompatible_ll
#   
#   N = sum(weights)
#   # I: number of incompatible genotypes
#   I = N - C
#   
#   
#   log_q_e = min(tmp$probs)
#   
#   incompatible_ll = 0.0
#   if( I > 0 ){
#     incompatible_ll = I * (log(1-alpha) + log_q_e)
#   }
#   
#   incompatible_ll = max(incompatible_ll, tmp2)
  
  list(ll=incompatible_ll + compatible_ll, incompatible_ll=incompatible_ll,  compatible_ll=compatible_ll, alpha=alpha)  
}
