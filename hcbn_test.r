rm(list=ls())
library(graph)
library(mccbn)
library(ggplot2)

# Set seed for reproducibility
set.seed(10)

################## INPUT OPTIONS ##################
p = 5                      # number of events
poset = random_poset(p)    # true poset
lambda_s = 1               # sampling rate
N = 100                    # number of observations / genotypes
eps = 0.05
##################################################

hamming_dist <- function(x, y) {
  dist = sum(x != y)
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

geno_prob_empirical <- function(N, poset, lambdas, lambda_s, genotype, eps) {
  # Compute P(genotype) empirically
  simGenotypes = sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas,
                                  eps=eps)
  return(sum(apply(simGenotypes$obs_events, 1, 
                   function(x, y) all(x == y), y=genotype)) / N)
}

# why samples are generated with noise? 
tdiff_empirical <- function(N, poset, lambdas, lambda_s, genotype, eps) {
  # Compute the time differences expirically for a given genotype
  simGenotypes = sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas,
                                  eps=eps)
  idx = which(apply(simGenotypes$obs_events, 1, function(x, y) all(x == y), 
                    y=genotype))
  
  return(apply(simGenotypes$T_events[idx, ], 2, mean))
}

dist_empirical <- function(N, poset, lambdas, lambda_s, genotype, eps = eps) {
  # Compute the average distance between genotype Y (subject to noise) and N
  # possible true genotypes X
  p = ncol(poset)
  simGenotypes = sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas, eps=eps)
  
  idx = which(apply(simGenotypes$obs_events, 1, function(x, y) all(x == y),  y=genotype))
  
  return(sum(apply(simGenotypes$hidden_genotypes[idx, ], 1, hamming_dist, y=genotype)) / length(idx))
}

prob_imp <- function(L, poset, lambdas, lambda_s, genotype, eps) {
  # Compute Pr(genotype) using Monte Carlo sampling
  simGenotypes = sample_genotypes(L, poset, sampling_param=lambda_s, lambdas=lambdas)
  p = ncol(poset)
  probs = apply(simGenotypes$obs_events, 1, function(x, y, e, p) { 
    d = sum(x != y)
    e^d * (1-e)^(p-d)
  }, y=genotype, e=eps, p=p)
  return(sum(probs) / L)
}

tdiff_imp <- function(L, poset, lambdas, lambda_s, genotype, eps) {
  # Compute expected time differences for a given genotype using Monte Carlo sampling
  p = ncol(poset)
  simGenotypes = sample_genotypes(L, poset, sampling_param=lambda_s, lambdas=lambdas)
  dist = apply(simGenotypes$obs_events, 1, hamming_dist, y=genotype)
  probs = eps^dist * (1-eps)^(p-dist)

  return(apply(simGenotypes$T_events, 2, function(x) sum(x * probs)) / sum(probs))
}

dist_imp <- function(L, poset, lambdas, lambda_s, genotype, eps) {
  # Compute expected hamming distance between a given observations (genotype) and the
  # underlying/true genotype using Monte Carlo sampling
  p = ncol(poset)
  simGenotypes = sample_genotypes(L, poset, sampling_param=lambda_s, lambdas=lambdas)
  dist = apply(simGenotypes$obs_events, 1, hamming_dist, y=genotype)
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
  eps = runif(1)
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
    for(i in 1:nrow(obs_events)) {
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
    
    llhood = N * sum(log(lambdas)) - sum(lambdas * t(expected_Tdiff)) + 
      p * log(eps) * sum(expected_dist) + p * log(1-eps) * sum(p - expected_dist)
    
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

###################################################################################
### MAIN PROGRAM
###################################################################################
# generate p random mutation rates uniformly distributed between lambda_s/3 to 3lambda_s.  
lambdas = runif(p, 1/3*lambda_s, 3*lambda_s)

# Simulate genotypes and sequencing times consistent with poset and mutation rates
# Sampling times are generated assuming they are exponentially distributed with rate
# lambda_s 
simulated_obs = sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas,
                                 eps=eps)

plot_poset(poset)

prob_empirical = numeric(N)
prob_sampling  = numeric(N)
for (i in 1:N) {
  genotype = simulated_obs$obs_events[i, ]
  prob_empirical[i] = geno_prob_empirical(N=100000, poset, lambdas, lambda_s, genotype, 
                                          eps = eps)
  prob_sampling[i] = prob_imp(L=100, poset, lambdas, lambda_s, genotype, eps = eps)
}

df = data.frame(x = prob_empirical, y = prob_sampling)
pl = ggplot(df, aes(x = x, y = y))
pl = pl + geom_point() + 
  geom_abline(intercept = 0, slope = 1, colour="blue") + 
  labs(x = "\nEmpirical probability", y = "Estimated probability\n" ) +
  theme_bw() +
  theme(text=element_text(size=14))

pl

time_diff_empirical = matrix(0, ncol=p, nrow=N)
time_diff_sampling  = matrix(0, ncol=p, nrow=N)
for (i in 1:N) {
  genotype = simulated_obs$obs_events[i, ]
  time_diff_empirical[i, ] = tdiff_empirical(N=100000, poset, lambdas, lambda_s, genotype, eps = eps)
  time_diff_sampling[i, ] = tdiff_imp(L=100, poset, lambdas, lambda_s, genotype, eps)
}

df = data.frame(x = time_diff_empirical[, 2], y = time_diff_sampling[, 2])
pl = ggplot(df, aes(x = x, y = y))
pl = pl + geom_point() + 
  geom_abline(intercept = 0, slope = 1, colour="blue") + 
  labs(x=expression('Empirical'~Z[2]), y=expression('Estimated'~Z[2]~'')) +
  theme_bw() +
  theme(text=element_text(size=14))

pl

X = possible_genotypes(p)
X_comp = apply(X, 1, is_compatible, poset=poset)
X_comp = X[X_comp, ]
mean(apply(X_comp, 1, hamming_dist, y=genotype))
dist_empirical(N=100000, poset, lambdas, lambda_s, genotype, eps)
dist_imp(L=100000, poset, lambdas, lambda_s, genotype, eps)

ret = MCMC_hcbn(poset, simulated_obs$obs_events)
abs(ret$lambdas - lambdas)/lambdas
abs(ret$avg_lambdas - lambdas)/lambdas
