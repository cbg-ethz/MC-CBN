

################ computing likelihood using importance sampling

cbn_density_ <- function(T, rates) {
  p = ncol(T)
  apply(T, 1, function(x) { sum(dexp(x, rates, log=TRUE))}  )
}

loglike_importance_sampling <- function(poset, lambda, obs_events, sampling_times, nrOfSamples, weights, lambda_s, sampling_times_available)  {
  
  topo_path = my.topological.sort(poset)-1
  probs = c()
  approx_loglike = 0.0
  
  for(i in 1:nrow(obs_events) ) {
    T = .Call("drawHiddenVarsSamples", nrOfSamples, obs_events[i,], sampling_times[i], lambda,  poset, topo_path, lambda_s, sampling_times_available)
    values = cbn_density_(T$T, lambda) - T$densities
    probs = c(probs, weights[i] * (logsumexp(values)$logR - log(nrOfSamples)))
    approx_loglike = approx_loglike + weights[i] * (logsumexp(values)$logR - log(nrOfSamples))
  }
  list(approx_loglike = approx_loglike , probs=probs)
}


genotype_prob_for_empty_poset <- function(lambda, genotype, time, lambda_s, sampling_time_available, log.p = FALSE) {
  p = length(lambda)
  
  if(sampling_time_available) {
    mu_probs = sapply(1:p, function(i) {pexp(time, lambda[i], lower.tail = genotype[i] == 1, log.p=TRUE) } )
  } else {
    mu_probs = sapply(1:p, function(i) {
      ratio = lambda[i]/(lambda[i] + lambda_s)
      ifelse(genotype[i] == 1, log(ratio), log(1-ratio))
      })
  }

  log_sum = sum(mu_probs)
  if(log.p == TRUE)
    return(log_sum)
  exp(log_sum)
}

all_genotype_prob_for_empty_poset <- function(lambda, obs_events, sampling_times, lambda_s, sampling_time_available, log.p=FALSE ) {
  
  probs = sapply(1:nrow(obs_events), function(i) { 
    genotype_prob_for_empty_poset(lambda, obs_events[i, ], sampling_times[i], lambda_s, sampling_time_available,  log.p = log.p)
  } )
  probs
}

#' Probability of a genotype given a poset and rate parameters, i.e. lambda
#' @export
genotype_probability_fast <- function(poset, lambda, genotype, time, nrOfSamples=100, lambda_s = 1.0, sampling_time_available=TRUE, topo_path=NULL) {
  
  if(is.null(topo_path) ) {
    topo_path = my.topological.sort(poset) -1
  }
  
  if(is_compatible(genotype, poset) == FALSE) {
    print(paste("The genotype (", paste(which(genotype==1), collapse=',') ,") is not compatible with the poset. Hence, the probability is zero") )
    return(0.0)
  }
  
  T = .Call("drawHiddenVarsSamples", nrOfSamples, genotype, time, lambda,  poset, topo_path, lambda_s, sampling_time_available)
  values = cbn_density_(T$T, lambda) - T$densities
  logp = (logsumexp(values)$logR - log(nrOfSamples))
  exp(logp)
}



########### likelihood


#' loglike_mixture_model
#' @export
loglike_mixture_model <- function(poset, lambda, obs_events, sampling_times, weights,  nrOfSamples, compatible_geno, 
                                  incomp_loglike, lambda_s, sampling_times_available){
  
  incompatible_ll = incomp_loglike$ll
  C = sum(weights[compatible_geno$compatible_indexes])
  alpha = incomp_loglike$alpha
  
  compatible_ll = 0.0
  if( C > 0 )
  {
    genotypes = obs_events[compatible_geno$compatible_indexes, , drop=FALSE ]
    
    sampling_times = sampling_times[compatible_geno$compatible_indexes]
    weights = weights[compatible_geno$compatible_indexes]
    
    tmp = loglike_importance_sampling(poset, lambda, genotypes, sampling_times, nrOfSamples, weights, lambda_s, sampling_times_available)
    compatible_ll = tmp$approx_loglike

    compatible_ll = compatible_ll + C * log(alpha)
  }
  
  list(ll=incompatible_ll + compatible_ll, incompatible_ll=incompatible_ll,  compatible_ll=compatible_ll, alpha=alpha)  
}
