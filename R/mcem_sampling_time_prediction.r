
################### Expected sampling time by importance sampling ######################

expected_sampling_time_importance_sampling <- function(genotype, poset, lambda, nrOfSamples ) {
  topo_path = my.topological.sort(poset) -1
  stimes = rexp(nrOfSamples, 1)
  
  T = list(T=matrix(0, 0, ncol(poset)), densities=c())
  for(stime in stimes) {
    t_ = .Call("drawHiddenVarsSamples", 1, genotype, stime, lambda,  poset, topo_path)  
    
    T$T  = rbind(T$T, t_$T)
    T$densities = c(T$densities, t_$densities)
  }
  
  values = cbn_density_(T$T, lambda) - T$densities
  #   logp = logsumexp(values)$logR
  
  #    weights = exp(cbn_density_(T$T, lambda) - T$densities)
  
  #   (sum(stimes * weights)/ nrOfSamples) /exp(logp) 
  #     sum(stimes * weights)/ sum(weights)
  
  logp = logsumexp(values)$logR
  exp(logsumexp(log(stimes) + values)$logR - logp)
}
