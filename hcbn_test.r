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
  simGenotypes = sample_genotypes(N, poset, sampling_param = lambda_s, lambdas=lambdas, eps=eps)
  sum(apply(simGenotypes$obs_events, 1, function(x) { all(x == genotype) } ))/N
}



genotype = simGenotypes$obs_events[5, ]
eps = 0.05

geno_prob_empirical(N=100000, poset, lambdas, lambda_s, genotype, eps = eps)


imp_prob <- function(L, poset, lambdas, lambda_s, genotype, eps) {
  simGenotypes = sample_genotypes(L, poset, sampling_param = lambda_s, lambdas=lambdas)
  p = ncol(poset)
  probs = apply(simGenotypes$obs_events, 1, function(x) { 
    d = sum(x != genotype)
    eps^d * (1-eps)^(p-d)
    } )
  sum(probs)/L
}

imp_prob(L=100, poset, lambdas, lambda_s, genotype, eps = eps)




