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
mccbn_path = "/Users/susanap/Documents/software/MC-CBN"
##################################################

source(file.path(mccbn_path, "hcbn_functions.r"))

empirical_vs_sampling <- function(events, L, rep=NULL, one_genotype=FALSE, 
                                  outname="") {
  
  if (one_genotype) {
    genotype = events
    N = rep
  } else {
    N = nrow(events)
  }
  
  prob_empirical = numeric(N)
  prob_sampling  = numeric(N)
  
  for (i in 1:N) {
    if (!one_genotype) {
      genotype = simulated_obs$obs_events[i, ]
    }
    prob_empirical[i] = geno_prob_empirical(N=100000, poset, lambdas, lambda_s, 
                                            genotype, eps=eps)
    prob_sampling[i] = prob_importance_sampling(genotype, L=L, poset, lambdas, 
                                lambda_s, eps=eps)
  }
  
  outname = file.path("~/Documents/hivX/CBN/hcbn_sampling/output/", 
                      paste("L", L, sep=""), 
                      paste("probability_Y_empirical_vs_sampling" , outname, 
                            ".pdf", sep=""))
  
  if (one_genotype) {
    df = data.frame(x = c(prob_empirical, prob_sampling), 
                    type = c(rep("empirical", N), rep("sampling", N)))
    pl = ggplot(df, aes(x = x, fill = type))
    pl = pl + geom_histogram(binwidth=0.01, alpha=0.5, position="identity") + 
      labs(x=expression(P(Y)), y="Frequency") + 
      theme_bw() + theme(text=element_text(size=14))
    
  } else {
    df = data.frame(x = prob_empirical, y = prob_sampling)
    pl = ggplot(df, aes(x = x, y = y))
    pl = pl + geom_point() + 
      geom_abline(intercept = 0, slope = 1, colour="blue") + 
      xlab(expression(P[empirical](Y))) + ylab(expression(widehat(P)(Y)~"\n")) +
      theme_bw() + theme(text=element_text(size=14))
    
  }
  ggsave(outname, pl, width=5, height=3)
  
  return(list("empirical" = prob_empirical, "sampling" = prob_sampling))
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
# NOTE: For p = 5, there are 32 possible genotypes. However, 16 are compatible with
#       poset

###################################################################################
### TEST 1
###################################################################################
# Check how well the sampling scheme approximates P(Y), Expected time differences and
# expected distance
###################### P(Y) ######################
L = 1000
probs = empirical_vs_sampling(simulated_obs$obs_events, L=L, 
                              outname=paste("_L", L, sep=""))

# One genotype
genotype = simulated_obs$obs_events[1, ] 
probs = empirical_vs_sampling(genotype, L=L, rep=N, one_genotype=TRUE, 
                              outname=paste("_g1_L", L, sep=""))

# WT
genotype = rep(0, p)
probs = empirical_vs_sampling(genotype, L=L, rep=N, one_genotype=TRUE, 
                              outname=paste("_WT_L", L, sep=""))

# Resistant type
genotype = rep(1, p)
probs = empirical_vs_sampling(genotype, L=L, rep=N, one_genotype=TRUE, 
                              outname=paste("_RT_L", L, sep=""))
# P(Y)
prob_empirical = numeric(N)
prob_sampling  = numeric(N)
for (i in 1:N) {
  genotype = simulated_obs$obs_events[i, ]
  prob_empirical[i] = geno_prob_empirical(N=100000, poset, lambdas, lambda_s, genotype, 
                                          eps=eps)
  prob_sampling[i] = prob_importance_sampling(genotype, L=100, poset, lambdas, lambda_s, eps=eps)
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
write(c(N, p+1), file.path(datadir, paste(filename, ".pat", sep="")))
write.table(cbind(rep(1, N), simulated_obs$obs_events), 
            file.path(datadir, paste(filename, ".pat", sep="")),
            row.names=FALSE, col.names=FALSE, append=TRUE)

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
write(c(N, p+1), file.path(datadir, paste(filename, ".pat", sep="")))
write.table(cbind(rep(1, N), simulated_obs$obs_events), 
            file.path(datadir, paste(filename, ".pat", sep="")),
            row.names=FALSE, col.names=FALSE, append=TRUE)
system(paste(hcbn_path, "h-cbn -f ", datadir, filename, " -w -v > ", datadir,
             filename, ".out.txt", sep=""))
lambdas_hcbn = read.csv(file.path(datadir, paste(filename, ".lambda", sep="")))
lambdas_hcbn = as.vector(t(lambdas_hcbn))
abs(lambdas_hcbn - lambdas)/lambdas

# sanity check (h-cbn report (N = 100)-302.451, (N = 1000) -3099.03)
obs_log_likelihood(simulated_obs$obs_events, poset, lambdas_hcbn, lambda_s, 
                   0.000036, L=10000) # N=100, -302.5706; N=1000, -3101.905
