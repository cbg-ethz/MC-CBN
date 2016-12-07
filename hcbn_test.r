rm(list=ls())
library(graph)
library(mccbn)

set.seed(10)


p = 5
poset = random_poset(p)
lambda_s = 1
N = 5
# generate p random mutation rates uniformly distributed between lambda_s/3 to 3lambda_s.  
lambdas = runif(p, 1/3*lambda_s, 3*lambda_s)
simGenotypes = sample_genotypes(N, poset, sampling_param = lambda_s, lambdas=lambdas)


plot_poset(poset)


geno_prob_empirical <- function(N, poset, lambdas, lambda_s, genotype, eps) {
  simGenotypes = sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas, eps=eps)
  sum(apply(simGenotypes$obs_events, 1, function(x) { all(x == genotype) } )) / N
}

tdiff_empirical <- function(N, poset, lambdas, lambda_s, genotype, eps) {
  simGenotypes = sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas, eps=eps)
  idx = which(apply(simGenotypes$obs_events, 1, function(x) { all(x == genotype) } ))
  
  return(apply(simGenotypes$T_events[idx, ], 2, mean))
}

genotype = simGenotypes$obs_events[1, ]
eps = 0.05

geno_prob_empirical(N=100000, poset, lambdas, lambda_s, genotype, eps = eps)
tdiff_empirical(N=100000, poset, lambdas, lambda_s, genotype, eps = eps)

imp_prob <- function(L, poset, lambdas, lambda_s, genotype, eps) {
  simGenotypes = sample_genotypes(L, poset, sampling_param = lambda_s, lambdas=lambdas)
  p = ncol(poset)
  probs = apply(simGenotypes$obs_events, 1, function(x) { 
    d = sum(x != genotype)
    eps^d * (1-eps)^(p-d)
    } )
  sum(probs)/L
}

tdiff_imp <- function(L, poset, lambdas, lambda_s, genotype, eps) {
  simGenotypes = sample_genotypes(L, poset, sampling_param = lambda_s, lambdas=lambdas)
  p = ncol(poset)
  probs = apply(simGenotypes$obs_events, 1, function(x, genotype,eps) { 
    d = sum(x != genotype)
    eps^d * (1-eps)^(p-d)
  }, genotype=genotype, eps=eps)
  return(apply(simGenotypes$T_events, 2, function(x) sum(x * probs)) / sum(probs))
  
}

imp_prob(L=100, poset, lambdas, lambda_s, genotype, eps = eps)
tdiff_imp(L=100, poset, lambdas, lambda_s, genotype, eps)



################

MCMC_hcbn <- function(poset, genotypes, sampling_times=NULL, max_iter=100,  zeta = 0.2, nrOfSamples = 5, verbose = TRUE, maxLambdaValue=10^6, lambda_s=1.0)  {
  p = ncol(poset) # number of mutations
  
  
  # initialize 
  new_eps
  new_lambda
  
  for(i in 1:maxIter) {
    
      # one  EM iteration
    # E step
    for(i in 1:nrow(genotypes)) {
      simGenotypes = sample_genotypes(L, .... new_lambda)
      hamming_distances = diff(X, Y) for all simGenotypes
      probs_Y_X = apply(simGenotypes$obs_events, 1 ..., new_eps)
      Z = simGenotypes$T_events[, mut]
      
      expectedD = sum(probs_Y_X * hamming_distances)/ sum(probs_Y_X)
      expectT = 0.0
      
      # rbind
    }
      
    # M step
    new_eps = mean(expectedD) / p
    
    new_lambda = 1 / apply(expectedT, 2, mean)
  }
  list(new_lambda, new_eps)
}

