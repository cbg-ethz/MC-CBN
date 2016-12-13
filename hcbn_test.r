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
N = 1000                   # number of observations / genotypes
eps = 0.05
hcbn_path = "/Users/susanap/Documents/software/ct-cbn-0.1.04b/"
datadir = "/Users/susanap/Documents/hivX/CBN/hcbn_sampling/testdata/"
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

prob_imp <- function(genotype, L, poset, lambdas, lambda_s, eps) {
  # Compute Pr(genotype) using Monte Carlo sampling
  # Generate L samples from poset with parameters 'lambdas' and 'lambda_s'. In 
  # particular epsilon is zero (default value)
  simGenotypes = sample_genotypes(L, poset, sampling_param=lambda_s, lambdas=lambdas)
  p = ncol(poset)
  d = apply(simGenotypes$hidden_genotypes, 1, hamming_dist, y=genotype)
  probs = eps^d * (1-eps)^(p-d)
  return(sum(probs) / L)
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

################ TEST1 ################
# P(Y)
prob_empirical = numeric(N)
prob_sampling  = numeric(N)
for (i in 1:N) {
  genotype = simulated_obs$obs_events[i, ]
  prob_empirical[i] = geno_prob_empirical(N=100000, poset, lambdas, lambda_s, genotype, 
                                          eps = eps)
  prob_sampling[i] = prob_imp(genotype, L=100, poset, lambdas, lambda_s, eps=eps)
}

df = data.frame(x = prob_empirical, y = prob_sampling)
pl = ggplot(df, aes(x = x, y = y))
pl = pl + geom_point() + 
  geom_abline(intercept = 0, slope = 1, colour="blue") + 
  labs(x = "\nEmpirical probability", y = "Estimated probability\n" ) +
  theme_bw() +
  theme(text=element_text(size=14))

pl

# Time differences 
time_diff_empirical = matrix(0, ncol=p, nrow=N)
time_diff_sampling  = matrix(0, ncol=p, nrow=N)
for (i in 1:N) {
  genotype = simulated_obs$obs_events[i, ]
  time_diff_empirical[i, ] = tdiff_empirical(N=100000, poset, lambdas, lambda_s,
                                             genotype, eps=eps)
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

# Hamming distance 
X = possible_genotypes(p)
X_comp = apply(X, 1, is_compatible, poset=poset)
X_comp = X[X_comp, ]
mean(apply(X_comp, 1, hamming_dist, y=genotype))

d_empirical = numeric(N)
d_sampling  = numeric(N)
for (i in 1:N) {
  genotype = simulated_obs$obs_events[i, ]
  d_empirical[i] = dist_empirical(N=100000, poset, lambdas, lambda_s, genotype,
                                  eps=eps)
  d_sampling[i] = dist_imp(L=100, poset, lambdas, lambda_s, genotype, eps)
}

df = data.frame(x = d_empirical, y = d_sampling)
pl = ggplot(df, aes(x = x, y = y))
pl = pl + geom_point() + 
  geom_abline(intercept = 0, slope = 1, colour="blue") + 
  labs(x = "\nEmpirical distance", y = "Estimated distance\n" ) +
  theme_bw() +
  theme(text=element_text(size=14))

pl

################ TEST2 ################
set.seed(10)
dist = rowSums(simulated_obs$obs_events != simulated_obs$hidden_genotypes)
llhood = complete_log_likelihood(lambdas, simulated_obs$T_events, dist, eps)
obs_log_likelihood(simulated_obs$obs_events, poset, lambdas, lambda_s, 
                   eps, L=10000) #-267.8936

# MC-CBN, error model: h-cbn
ret = MCMC_hcbn(poset, simulated_obs$obs_events)
abs(ret$lambdas - lambdas)/lambdas
abs(ret$avg_lambdas - lambdas)/lambdas
abs(ret$eps - eps)
abs(ret$llhood - llhood)
obs_log_likelihood(simulated_obs$obs_events, poset, ret$avg_lambdas, lambda_s, 
                   ret$avg_eps, L=10000) ## N=100, L=100, -268.4495

# MC-CBN, error model: mixture model
compatible_idx = apply(simulated_obs$obs_events, 1, is_compatible, poset=poset)
ret_mixture = estimate_mutation_rates(poset, 
                                      simulated_obs$obs_events[compatible_idx, ])
abs(ret_mixture$par - lambdas)/lambdas
ret_mixture$ll ## N=100, L=5, -432.9113

# hcbn
filename = paste("simulated_obs_n", N, "_p", p, sep="")
write.csv(cbind(rep(1, N), simulated_obs$obs_events), 
          file.path(datadir, paste(filename, ".pat", sep="")),
          row.names=FALSE)
system(paste(hcbn_path, "h-cbn -f", datadir, filename, " -w -v > ", datadir,
             filename, ".out.txt", sep=""))
lambdas_hcbn = read.csv(file.path(datadir, paste(filename, ".lambda", sep="")))
lambdas_hcbn = as.vector(t(lambdas_hcbn))
abs(lambdas_hcbn - lambdas)/lambdas

# sanity check (h-cbn report -267.308)
obs_log_likelihood(simulated_obs$obs_events, poset, lambdas_hcbn, lambda_s, 
                   0.034023, L=10000) # N=100, -268.256

################ TEST3 ################
set.seed(10)
# empty poset
poset = matrix(0, p, p)

simulated_obs = sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas,
                                 eps=eps)

obs_log_likelihood(simulated_obs$obs_events, poset, lambdas, lambda_s, eps, 
                   exact=TRUE) # N=100, -340.4663; N=1000, -3374.1

# MC-CBN, error model: h-cbn
ret = MCMC_hcbn(poset, simulated_obs$obs_events)
obs_log_likelihood(simulated_obs$obs_events, poset, ret$avg_lambdas, lambda_s, 
                   ret$avg_eps, L=10000) # N=100, -305.5238; N=1000, -3109.656

filename = paste("simulated_obs_n", N, "_p", p, "_empty", sep="")
write.csv(cbind(rep(1, N), simulated_obs$obs_events), 
          file.path(datadir, paste(filename, ".pat", sep="")),
          row.names=FALSE)
system(paste(hcbn_path, "h-cbn -f ", datadir, filename, " -w -v > ", datadir,
             filename, ".out.txt", sep=""))
lambdas_hcbn = read.csv(file.path(datadir, paste(filename, ".lambda", sep="")))
lambdas_hcbn = as.vector(t(lambdas_hcbn))
abs(lambdas_hcbn - lambdas)/lambdas

# sanity check (h-cbn report (N = 100)-302.451, (N = 1000) -3099.03)
obs_log_likelihood(simulated_obs$obs_events, poset, lambdas_hcbn, lambda_s, 
                   0.000036, L=10000) # N=100, -302.5706; N=1000, -3101.905
