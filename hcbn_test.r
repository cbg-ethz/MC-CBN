rm(list=ls())
library(graph)
library(mccbn)
library(ggplot2)

# Set seed for reproducibility
set.seed(10)

############################### INPUT OPTIONS ################################
p = 5                      # number of events
poset = random_poset(p)    # true poset
lambda_s = 1               # sampling rate
N = 100 #1000              # number of observations / genotypes
eps = 0.05
hcbn_path = "/Users/susanap/Documents/software/ct-cbn-0.1.04b/"
datadir = "/Users/susanap/Documents/hivX/CBN/hcbn_sampling/testdata/"
outdir = "/Users/susanap/Documents/hivX/CBN/hcbn_sampling/output/" 
mccbn_path = "/Users/susanap/Documents/software/MC-CBN"
###############################################################################

source(file.path(mccbn_path, "hcbn_functions.r"))

###############################################################################
### MAIN PROGRAM
###############################################################################
# generate p random mutation rates uniformly distributed between lambda_s/3 to 
# 3 lambda_s.  
lambdas = runif(p, 1/3*lambda_s, 3*lambda_s)

# Simulate genotypes and sequencing times consistent with poset and mutation rates
# Sampling times are generated assuming they are exponentially distributed with 
# rate lambda_s 
simulated_obs = sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas,
                                 eps=eps)

plot_poset(poset)
# NOTE: For p = 5, there are 32 possible genotypes. However, 16 are compatible 
#       with poset

###############################################################################
### TEST 1
###############################################################################
# Check how well the sampling scheme approximates P(Y), Expected time 
# differences and expected distance
#################################### P(Y) #####################################
L = 1000
binwidth = 0.01
# binwidth: L=10 0.1; L=100 ; L=1000 0.01
probs = prob_empirical_vs_sampling(simulated_obs$obs_events, L=L,  
                                   sampling='naive', outdir=outdir, 
                                   outname=paste("_N", N, "_L", L, sep=""))
save(probs, file=file.path(outdir, paste("L", L, "_naive", sep=""), 
                           paste("probability_Y_N", N, "_L", L, ".RData",
                                 sep="")))

# One genotype
genotype = simulated_obs$obs_events[1, ] 
probs = prob_empirical_vs_sampling(genotype, L=L, rep=N, one_genotype=TRUE, 
                                   sampling='naive', outdir=outdir, 
                                   outname=paste("_g1_rep", N, "_L", L, sep=""), 
                                   binwidth=binwidth)
save(probs, file=file.path(outdir, paste("L", L, "_naive", sep=""), 
                           paste("probability_Y_g1_rep", N, "_L", L, ".RData",
                                 sep="")))

# WT
genotype = rep(0, p)
probs = prob_empirical_vs_sampling(genotype, L=L, rep=N, one_genotype=TRUE, 
                                   sampling='naive', outdir=outdir,
                                   outname=paste("_WT_rep", N, "_L", L, sep=""),
                                   binwidth=binwidth)
save(probs, file=file.path(outdir, paste("L", L, "_naive", sep=""), 
                           paste("probability_Y_WT_rep", N, "_L", L, ".RData",
                                 sep="")))

# Resistant type
genotype = rep(1, p)
probs = prob_empirical_vs_sampling(genotype, L=L, rep=N, one_genotype=TRUE, 
                                   sampling='naive', outdir=outdir,
                                   outname=paste("_RT_rep", N, "_L", L, sep=""),
                                   binwidth=binwidth)
save(probs, file=file.path(outdir, paste("L", L, "_naive", sep=""), 
                           paste("probability_Y_RT_rep", N, "_L", L, ".RData",
                                 sep="")))

########################## Expected time differences ##########################
set.seed(10)
time_diff = tdiff_empirical_vs_sampling(simulated_obs$obs_events, L=L,  
                                        sampling='naive', outdir=outdir, 
                                        outname=paste("_N", N, "_L", L, sep=""))
save(time_diff, file=file.path(outdir, paste("L", L, "_naive", sep=""), 
                               paste("time_diff_N", N, "_L", L, ".RData", 
                                     sep="")))

############################## Expected distance ############################## 
set.seed(10)
dist = dist_empirical_vs_sampling(simulated_obs$obs_events, L=L, 
                                  sampling='naive', outdir=outdir, 
                                  outname=paste("_N", N, "_L", L, sep=""))
save(dist, file=file.path(outdir, paste("L", L, "_naive", sep=""), 
                          paste("hamming_dist_N", N, "_L", L, ".RData", sep="")))

X = possible_genotypes(p)
X_comp = apply(X, 1, is_compatible, poset=poset)
X_comp = X[X_comp, ]
mean(apply(X_comp, 1, hamming_dist, y=genotype))


###############################################################################
### TEST 2
###############################################################################
# Check that observed log-likelihood is similar among truth (approximate), 
# MC-CBN using H-CBN error model and mixture model, and H-CBN
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

###############################################################################
### TEST 3
###############################################################################
# Comparison between truth (exact) and MC-CBN using H-CBN error model. For the
# empty poset, observed log-likelihood can be computed exactly
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


###############################################################################
### TEST 3
###############################################################################
# Effective sample size for different poset samples
library(doMC)
thrds = 4
registerDoMC(thrds)

set.seed(10)
p = 2^seq(2, 5, 1)
N = rep(50, 4) #sapply(50*p, min, 1000)
L = 1000

effective_sample_size = foreach (k = 1:length(p)) %dopar% {
  poset = random_poset(p[k])
  lambdas = runif(p[k], 1/3*lambda_s, 3*lambda_s)
  simulated_obs = sample_genotypes(N[k], poset, sampling_param=lambda_s, 
                                   lambdas=lambdas, eps=eps)
  Ne = numeric(N[k])
  for (i in 1:N[k]) {
    w = importance_weight(simulated_obs$obs_events[i, ], L=L, poset, lambdas,
                          lambda_s, eps, sampling='naive')
    Ne[i] = sum(w$w)^2 / sum(w$w^2)
  }
  return(Ne)
}

df = data.frame(x=rep(p, N), y=unlist(effective_sample_size))

pl = ggplot(df, aes(x=factor(x), y=y)) + 
  geom_boxplot(varwidth=TRUE, fill='cornflowerblue') + 
  labs(x="\n poset size", y="effective sample \n") + 
  theme_bw() + theme(text=element_text(size=14))

ggsave(file.path(outdir, paste("L", L, "_naive", sep=""), 
                 "effective_sample_size.pdf"),
       pl, width=3.5, height=2.5)
