
########## comparison of time differences

time_differences_empirical <- function(N, poset, lambdas, lambda_s, genotype, eps) {
  simGenotypes = sample_genotypes(N, poset, sampling_param = lambda_s, lambdas=lambdas, eps=eps)
  
  indexes = which(apply(simGenotypes$obs_events, 1, function(x) { all(x == genotype) } ) )
  
  apply(simGenotypes$T_events[indexes, ], 2, mean )
}


time_differences_empirical(N=100000, poset, lambdas, lambda_s, genotype, eps = eps)

plot_poset(poset)
1/lambdas 

