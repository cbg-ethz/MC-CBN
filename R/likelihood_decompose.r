
.log_observed_geno_probs_decomp <- function(poset, lambda, obs_events, dtimes, with_eps=FALSE, eps=NA) {
  graphobj <- graph.adjacency(poset, mode="undirected")
  membership = clusters(graphobj)$membership
  
  lprobs = rep(0, nrow(obs_events) )
  for(i in 1:max(membership)) {
    indexes = which(membership == i)
    G = geno_indexes = NA
    if(length(indexes) > 1) {
      G = compatible_genotypes_with_matrix(poset[indexes, indexes])
      geno_indexes = compute_geno_indexes(obs_events[, indexes], poset[indexes, indexes])  
    } 
    if(with_eps==TRUE) {
      G = compatible_genotypes_with_matrix(poset[indexes, indexes])
      lprobs = lprobs + log_observed_geno_probs_with_eps(lambda[indexes], eps, obs_events[, indexes], 
                                                dtimes, G)
    } else{ 
      lprobs = lprobs + log_observed_geno_probs(lambda[indexes], obs_events[, indexes], 
                                                dtimes, G, geno_indexes)
    }
    
  }
  lprobs
}



.likelihood_decomp <- function(poset, lambda, obs_events, dtimes) {
  sum(.log_observed_geno_probs_decomp(poset, lambda, obs_events, dtimes))  
}




loglike <- function(poset, lambda, obs_events, times, D, with_eps=FALSE, eps=NA)  {
  dtimes = discretize_times(times, D)
  
  sum(.log_observed_geno_probs_decomp(poset, lambda, obs_events, dtimes, with_eps, eps) )
}