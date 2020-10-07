rm(list=ls())
library(graph)
library(mccbn)
library(ggplot2)

# Script for performing several tests/checks
# 1. Check how well the sampling scheme approximates P(Y), expected time 
#    differences and expected distance (Emprirical vs sampling)
# 2. Check that observed log-likelihood is similar among (i) truth (approx.), 
#    (ii) MC-CBN using H-CBN error model and (iii) MC-CBN using mixture error
#    model, and (iv) H-CBN.
# 3. Comparison between truth (exact), MC-CBN using H-CBN error model, MC-CBN
#    using mixture error model and H-CBN for an empty poset. For empty posets, 
#    observed log-likelihood can be computed exactly
# 4. Effective sample size for a particular observation, as well as multiple
#    observations
# 5. Effective sample size for different poset sizes
# 6. Effective sample size as function of (i) compatible/incompatible status, 
#    (ii) number of mutations (1's) per observation, (iii) number of mutations 
#    that could be added or (iv) removed, and (iv) hamming distance 

# Set seed for reproducibility
set.seed(10)

############################### INPUT OPTIONS ################################
p = 5                      # number of events
poset = random_poset(p)    # true poset
lambda_s = 1               # sampling rate
n = 100                    # number of observations / genotypes
L = 1000                   # number of samples
eps = 0.05                 # true epsilon
sampling = 'naive'         # options: 'naive', 'add-remove'
perturb_prob = 0           # option used if sampling is set to 'add-remove'
run_hcbn = FALSE           # indicate whether or not H-CBN should be executed
save_output = FALSE        # indicate whether or not save output

# The following options are NOT MANDATORY to specify/modify
# Modify the following options, if you want to perform comparisons to h-cbn 
# (if so, option 'run_hcbn' should be set to TRUE)
hcbn_path = "/Users/susanap/Documents/software/ct-cbn-0.1.04b/"
datadir = "/Users/susanap/Documents/hivX/CBN/hcbn_sampling/testdata/"

# Specify the directory where output files are to be saved. It is needed, if 
# you want to save output files (If so, option 'save_output' should be set to
# TRUE). If path doesn't exit, it will set to the working directory
outdir = "/Users/susanap/Documents/hivX/CBN/hcbn_sampling/output"

# Specify where to find 'hcbn_functions.r' (directory where MC-CBN library is
# saved). If path doesn't exit, it will set to the working directory
mccbn_path = "/Users/susanap/Documents/software/MC-CBN"

###############################################################################

### Check that required paths are given and exist
if (!dir.exists(mccbn_path)) {
  cat("Specified directory doesn't exist. Setting 'mccbn_path' to working" ,
      " directory, \'", getwd(), "\'\n", sep="")
  mccbn_path = getwd()
} 
if (!save_output) {
  outdir = NULL
} else if (!dir.exists(outdir) & save_output) {
  cat("Specified directory doesn't exist. Setting 'outdir' to working" ,
      " directory, \'", getwd(), "\'\n", sep="")
  outdir = getwd()
} 
if (!dir.exists(hcbn_path) & run_hcbn) {
  stop("Path to H-CBN executables should be specified")
}
if (!dir.exists(datadir) & run_hcbn) {
  cat("Specified directory doesn't exist. Setting 'datadir' to working" ,
      " directory, \'", getwd(), "\'\n", sep="")
  datadir = getwd()
}

### Load functions 
if (file.exists(file.path(mccbn_path, "hcbn_functions.r"))) {
  cat("Loading functions from \'" , mccbn_path, "\'\n", sep="")
  source(file.path(mccbn_path, "hcbn_functions.r"))
} else {
  stop("Script 'hcbn_functions.r' cannot be found in \'", mccbn_path, "\'")
}

###############################################################################
### MAIN PROGRAM
###############################################################################
# generate p random mutation rates uniformly distributed between lambda_s/3 to 
# 3 lambda_s.  
lambdas = runif(p, 1/3*lambda_s, 3*lambda_s)

# Simulate genotypes and sequencing times consistent with poset and mutation
# rates. Sampling times are generated assuming they are exponentially
# distributed with rate lambda_s
simulated_obs =
  sample_genotypes(n, poset, sampling_param=lambda_s, lambdas=lambdas, eps=eps)

plot_poset(poset)
# NOTE: For p = 5, there are 32 possible genotypes. However, 16 are compatible 
#       with poset

###############################################################################
### TEST 1
###############################################################################
# Check how well the sampling scheme approximates P(Y), expected time 
# differences and expected distance
if (L <= 100) {
  binwidth = 0.1
} else if (L <= 1000) {
  binwidth = 0.01
}

if (save_output) {
  out_dir = file.path(outdir, paste("L", L, "_", sampling, sep=""))
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
}

#################################### P(Y) #####################################
probs = prob_empirical_vs_sampling(simulated_obs$obs_events, L=L,  
                                   sampling=sampling, outdir=outdir, 
                                   outname=paste("_N", n, "_L", L, sep=""))
if (save_output) {
  save(probs, file=file.path(out_dir, 
                             paste("probability_Y_N", n, "_L", L, ".RData", 
                                   sep="")))
}

# One genotype
genotype = simulated_obs$obs_events[1, ] 
probs = prob_empirical_vs_sampling(genotype, L=L, rep=n, one_genotype=TRUE, 
                                   sampling=sampling, outdir=outdir, 
                                   outname=paste("_g1_rep", n, "_L", L, sep=""), 
                                   binwidth=binwidth)
if (save_output) {
  save(probs, file=file.path(out_dir, 
                             paste("probability_Y_g1_rep", n, "_L", L, ".RData", 
                                   sep="")))
}

# Wild type
genotype = rep(0, p)
probs = prob_empirical_vs_sampling(genotype, L=L, rep=n, one_genotype=TRUE, 
                                   sampling=sampling, outdir=outdir,
                                   outname=paste("_WT_rep", n, "_L", L, sep=""),
                                   binwidth=binwidth)
if (save_output) {
  save(probs, file=file.path(out_dir, 
                             paste("probability_Y_WT_rep", n, "_L", L, ".RData",
                                   sep="")))
}

# Resistant type
genotype = rep(1, p)
probs = prob_empirical_vs_sampling(genotype, L=L, rep=n, one_genotype=TRUE, 
                                   sampling=sampling, outdir=outdir,
                                   outname=paste("_RT_rep", n, "_L", L, sep=""),
                                   binwidth=binwidth)
if (save_output) {
  save(probs, file=file.path(out_dir,
                             paste("probability_Y_RT_rep", n, "_L", L, ".RData",
                                   sep="")))
}


########################## Expected time differences ##########################
set.seed(10)
time_diff = tdiff_empirical_vs_sampling(simulated_obs$obs_events, L=L,  
                                        sampling=sampling, outdir=outdir, 
                                        outname=paste("_N", n, "_L", L, sep=""))
if (save_output) {
  save(time_diff, file=file.path(out_dir, 
                                 paste("time_diff_N", n, "_L", L, ".RData", 
                                       sep="")))
}


############################## Expected distance ############################## 
set.seed(10)
dist = dist_empirical_vs_sampling(simulated_obs$obs_events, L=L, 
                                  sampling=sampling, outdir=outdir, 
                                  outname=paste("_N", n, "_L", L, sep=""))

if (save_output) {
  save(dist, file=file.path(out_dir, 
                            paste("hamming_dist_N", n, "_L", L, ".RData", sep="")))
}

X = possible_genotypes(p)
X_comp = apply(X, 1, is_compatible, poset=poset)
X_comp = X[X_comp, ]
mean(apply(X_comp, 1, hamming_dist, y=genotype))

###############################################################################
### TEST 2
###############################################################################
# Check that observed log-likelihood is similar among (i) truth (approximate), 
# (ii) MC-CBN using H-CBN error model and (iii) MC-CBN using mixture error model,
# and (iv) H-CBN

dist = rowSums(simulated_obs$obs_events != simulated_obs$hidden_genotypes)
## complete-data log-likelihood
llhood = complete_log_likelihood(lambdas, simulated_obs$T_events, dist, eps) 

## observed log-likelihood (true poset and params - lambdas, lambda_s and eps)
obs_log_likelihood(simulated_obs$obs_events, poset, lambdas, lambda_s, eps,
                   L=10000, sampling=sampling, seed=10) 

# MC-CBN, error model: h-cbn
# by default, L = 100
ret = MCEM_hcbn(poset, simulated_obs$obs_events, sampling_times=NULL, 
                lambda_s=lambda_s, sampling=sampling, parallel=FALSE, seed=10)
ret$llhood
ret$avg_llhood

obs_log_likelihood(simulated_obs$obs_events, poset, ret$avg_lambdas, lambda_s, 
                   ret$avg_eps, L=10000, sampling=sampling, seed=10)

# MC-CBN, error model: mixture model
set.seed(10)
compatible_obs = compatible_genotypes(simulated_obs$obs_events, poset)
ret_mixture = 
  estimate_mutation_rates(poset,
                          simulated_obs$obs_events[compatible_obs$compatible_indexes, ],
                          nrOfSamples=100)
ret_mixture$ll 

geno_prob_noise = 
  genotype_probs_empty_poset(
    simulated_obs$obs_events[-compatible_obs$compatible_indexes, ], 
    sampling_times=rep(0, (1 - compatible_obs$fraction) * n), 
    weights=rep(1, (1 - compatible_obs$fraction) * n), max_iter=100, zeta=0.2,
    nrOfSamplesForEStep=100, verbose=FALSE, maxLambdaValue=1e6, 
    lambda_s=lambda_s, sampling_times_available=FALSE)

llhood_incompatible = 
  incompatible_loglike(poset, simulated_obs$obs_events, sampling_times=rep(0, n), 
                       weights=rep(1, n), compatible_geno=compatible_obs, 
                       noise_model="empty", geno_prob_noise=geno_prob_noise) 
llhood_incompatible$ll
llhood_incompatible$alpha

loglike_mixture_model(poset, ret_mixture$par, simulated_obs$obs_events, 
                      sampling_times=rep(0, n), weights=rep(1, n), 
                      nrOfSamples=100, compatible_geno=compatible_obs, 
                      incomp_loglike=llhood_incompatible, lambda_s=lambda_s, 
                      sampling_times_available=FALSE) 

# hcbn
if (run_hcbn) {
  filename = paste("simulated_obs_n", n, "_p", p, sep="")
  if (!file.exists(file.path(datadir, paste(filename, ".pat", sep="")))) {
    write(c(n, p+1), file.path(datadir, paste(filename, ".pat", sep="")))
    write.table(cbind(rep(1, n), simulated_obs$obs_events), 
                file.path(datadir, paste(filename, ".pat", sep="")),
                row.names=FALSE, col.names=FALSE, append=TRUE)
  }
  if (!file.exists(file.path(datadir, paste(filename, ".poset", sep="")))) {
    write_poset(poset, filename, datadir)
  }
  if (!dir.exists(file.path(datadir,filename))) {
    dir.create(file.path(datadir,filename))
  }
  system(paste(hcbn_path, "h-cbn -f", datadir, filename, " -w > ", datadir,
               filename, ".out.txt", sep=""))
  params_hcbn = read.csv(file.path(datadir, paste(filename, ".out.txt", 
                                                  sep="")), sep="\t")
  lambdas_hcbn = as.numeric(params_hcbn[1, 6:ncol(params_hcbn)])
  eps_hcbn = params_hcbn$Eps
  
  params_hcbn$Loglik
  # sanity check (in order to compare with what hcbn reports)
  obs_log_likelihood(simulated_obs$obs_events, poset, lambdas_hcbn, 
                     params_hcbn$lambda_s, eps_hcbn, L=10000, 
                     sampling=sampling, seed=10) 
}

## COMPARISONS
# MC-CBN, hidden layer
abs(ret$lambdas - lambdas)/lambdas
median(abs(ret$lambdas - lambdas))/median(lambdas)
abs(ret$avg_lambdas - lambdas)/lambdas
median(abs(ret$avg_lambdas - lambdas))/median(lambdas)
abs(ret$eps - eps)
abs(ret$avg_eps - eps)
abs(ret$llhood - llhood)
abs(ret$avg_llhood - llhood)

# MC-CBN, mixture error model
abs(ret_mixture$par - lambdas)/lambdas
median(abs(ret_mixture$par - lambdas))/median(lambdas)

# H-CBN
if (run_hcbn) {
  abs(lambdas_hcbn - lambdas)/lambdas
  median(abs(lambdas_hcbn - lambdas))/median(lambdas)
  abs(eps_hcbn - eps)
}

###############################################################################
### TEST 3
###############################################################################
# Comparison between truth (exact) and MC-CBN using H-CBN error model. For the
# empty poset, observed log-likelihood can be computed exactly
set.seed(10)
# empty poset
poset_empty = matrix(0, p, p)

simulated_obs_empty = sample_genotypes(n, poset_empty, sampling_param=lambda_s, 
                                       lambdas=lambdas, eps=eps)

dist = rowSums(simulated_obs_empty$obs_events != simulated_obs_empty$hidden_genotypes)
## complete-data log-likelihood
llhood = complete_log_likelihood(lambdas, simulated_obs_empty$T_events, dist,
                                 eps) 

## observed log-likelihood (true poset and params - lambdas, lambda_s and eps)
obs_log_likelihood(simulated_obs_empty$obs_events, poset_empty, lambdas, 
                   lambda_s, eps, exact=TRUE) #n=1000, -3374.1
obs_log_likelihood(simulated_obs_empty$obs_events, poset_empty, lambdas, 
                   lambda_s, eps, L=10000, sampling=sampling, seed=10) 

# MC-CBN, error model: h-cbn
# by default, L = 100
ret = MCEM_hcbn(poset_empty, simulated_obs_empty$obs_events, sampling_times=NULL, 
                lambda_s=lambda_s, sampling=sampling, parallel=FALSE, seed=10)
ret$llhood
ret$avg_llhood

obs_log_likelihood(simulated_obs_empty$obs_events, poset_empty, ret$avg_lambdas,
                   lambda_s, ret$avg_eps, L=10000, sampling=sampling, seed=10)

ret$llhood
ret$avg_llhood

obs_log_likelihood(simulated_obs_empty$obs_events, poset_empty, ret$avg_lambdas,
                   lambda_s, ret$avg_eps, L=10000) # n=1000, -3109.656

# MC-CBN, error model: mixture model (all observations are compatible with
# empty poset)
set.seed(10)
compatible_obs = compatible_genotypes(simulated_obs_empty$obs_events, 
                                      poset_empty)
stopifnot(compatible_obs$fraction == 1)
ret_mixture = 
  estimate_mutation_rates(poset_empty, simulated_obs_empty$obs_events, 
                          nrOfSamples=100)
ret_mixture$ll 

llhood_incompatible = 
  incompatible_loglike(poset_empty, simulated_obs_empty$obs_events, 
                       sampling_times=rep(0, n), weights=rep(1, n), 
                       compatible_geno=compatible_obs, noise_model="empty",
                       geno_prob_noise=0) 
llhood_incompatible$ll
llhood_incompatible$alpha

loglike_mixture_model(poset_empty, ret_mixture$par, 
                      simulated_obs_empty$obs_events, sampling_times=rep(0, n),
                      weights=rep(1, n), nrOfSamples=100, 
                      compatible_geno=compatible_obs, 
                      incomp_loglike=llhood_incompatible, lambda_s=lambda_s, 
                      sampling_times_available=FALSE) 

if (run_hcbn) {
  filename = paste("simulated_obs_n", n, "_p", p, "_empty", sep="")
  if (!file.exists(file.path(datadir, paste(filename, ".pat", sep="")))) {
    write(c(n, p+1), file.path(datadir, paste(filename, ".pat", sep="")))
    write.table(cbind(rep(1, n), simulated_obs_empty$obs_events), 
                file.path(datadir, paste(filename, ".pat", sep="")),
                row.names=FALSE, col.names=FALSE, append=TRUE)
  }
  if (!file.exists(file.path(datadir, paste(filename, ".poset", sep="")))) {
    write_poset(poset_empty, filename, datadir)
  }
  if (!dir.exists(file.path(datadir,filename))) {
    dir.create(file.path(datadir,filename))
  }
  system(paste(hcbn_path, "h-cbn -f ", datadir, filename, " -w > ", datadir,
               filename, ".out.txt", sep=""))
  params_hcbn = read.csv(file.path(datadir, paste(filename, ".out.txt", 
                                                  sep="")), sep="\t")
  lambdas_hcbn = as.numeric(params_hcbn[1, 6:ncol(params_hcbn)])
  eps_hcbn = params_hcbn$Eps
  
  params_hcbn$Loglik
  # sanity check (in order to compare with what hcbn reports, e.g. for N = 1000, -3099.03)
  obs_log_likelihood(simulated_obs_empty$obs_events, poset_empty, lambdas_hcbn, 
                     params_hcbn$lambda_s, eps_hcbn, L=10000, 
                     sampling=sampling, seed=10)
  # N=1000, -3101.905
}

## COMPARISONS
# MC-CBN, hidden layer
abs(ret$lambdas - lambdas)/lambdas
median(abs(ret$lambdas - lambdas))/median(lambdas)
abs(ret$avg_lambdas - lambdas)/lambdas
median(abs(ret$avg_lambdas - lambdas))/median(lambdas)
abs(ret$eps - eps)
abs(ret$avg_eps - eps)
abs(ret$llhood - llhood)
abs(ret$avg_llhood - llhood)

# MC-CBN, mixture error model
abs(ret_mixture$par - lambdas)/lambdas
median(abs(ret_mixture$par - lambdas))/median(lambdas)

if (run_hcbn) {
  abs(lambdas_hcbn - lambdas)/lambdas
  median(abs(lambdas_hcbn - lambdas))/median(lambdas)
  abs(eps_hcbn - eps)
}

###############################################################################
### TEST 4 
###############################################################################
# EXPLORATORY: Effective sample size for particular observations and on average 
# for 100 obervations
sampling = 'add-remove'
idx = 2
genotype = simulated_obs$obs_events[idx, ]
t0 = Sys.time()
w = importance_weight(genotype=genotype, L=L, poset=poset, lambdas=lambdas,
                      lambda_s=lambda_s, eps=eps, sampling_time=NULL, 
                      perturb_prob=perturb_prob, sampling=sampling, seed=10)
Ne = sum(w$w)^2 / sum(w$w^2)
print(difftime(Sys.time(), t0))
print(Ne)

# contrasting with known sampling time
t0 = Sys.time()
w = importance_weight(genotype=genotype, L=L, poset=poset,
                      lambdas=lambdas, lambda_s=lambda_s, eps=eps,
                      sampling_time=simulated_obs$T_sampling[idx], 
                      perturb_prob=perturb_prob, sampling=sampling, seed=10)
Ne = sum(w$w)^2 / sum(w$w^2)
print(difftime(Sys.time(), t0))
print(Ne)

# Over more observations
t0 = Sys.time()
k = min(n, 100)
Ne = numeric(k)
for (i in 1:k) {
  w = importance_weight(genotype=simulated_obs$obs_events[i, ], L=L,
                        poset=poset, lambdas=lambdas, lambda_s=lambda_s,
                        eps=eps, perturb_prob=perturb_prob, sampling_time=NULL,
                        sampling=sampling, seed=10)
  Ne[i] = sum(w$w)^2 / sum(w$w^2)
}
cat(difftime(Sys.time(), t0, units='mins'), mean(Ne*100/L), sd(Ne*100/L), 
    max(Ne*100/L), min(Ne*100/L), sep=" ")

# contrasting with known sampling time
t0 = Sys.time()
Ne = numeric(100)
for (i in 1:k) {
  w = importance_weight(genotype=simulated_obs$obs_events[i, ], L=L,
                        poset=poset, lambdas=lambdas, lambda_s=lambda_s,
                        eps=eps, perturb_prob=perturb_prob, 
                        sampling_time=simulated_obs$T_sampling[i], 
                        sampling=sampling, seed=10)
  Ne[i] = sum(w$w)^2 / sum(w$w^2)
}
cat(difftime(Sys.time(), t0, units='mins'), mean(Ne*100/L), sd(Ne*100/L), 
    max(Ne*100/L), min(Ne*100/L), sep=" ")

###############################################################################
### TEST 5
###############################################################################
# Effective sample size for different poset sizes
library(doMC)
thrds = 4
registerDoMC(thrds)

set.seed(10)
p = 2^seq(2, 5, 1)
N = sapply(50*p, min, 1000)
L = 1000

t0 <- Sys.time()
effective_sample_size = foreach (k = 1:length(p)) %dopar% {
  poset = random_poset(p[k])
  lambdas = runif(p[k], 1/3*lambda_s, 3*lambda_s)
  simulated_obs = sample_genotypes(N[k], poset, sampling_param=lambda_s, 
                                   lambdas=lambdas, eps=eps)
  Ne = numeric(N[k])
  for (i in 1:N[k]) {
    w = importance_weight(simulated_obs$obs_events[i, ], L=L, poset, lambdas,
                          lambda_s, eps, perturb_prob=perturb_prob, 
                          sampling_time=NULL, sampling=sampling, seed=10)
    Ne[i] = sum(w$w)^2 / sum(w$w^2)
  }
  return(Ne)
}
runtime = as.numeric(difftime(Sys.time(), t0, units='mins'))

df = data.frame(x=rep(p, N), y=unlist(effective_sample_size) * 100 /L)

pl = ggplot(df, aes(x=factor(x), y=y)) + 
  geom_boxplot(varwidth=TRUE, fill='cornflowerblue') + 
  labs(x="\n poset size", y="effective sample [%] \n") + 
  theme_bw() + theme(text=element_text(size=14))

if (save_output) {
  ggsave(file.path(outdir, paste("L", L, "_", sampling, sep=""), 
                   "effective_sample_size.pdf"),
         pl, width=3.5, height=2.5)
}

############## Alternative: parallel execution of the sampling scheme
effective_sample_size = list()
t0 <- Sys.time()
for (k in 1:length(p)) {
  poset = random_poset(p[k])
  lambdas = runif(p[k], 1/3*lambda_s, 3*lambda_s)
  simulated_obs = sample_genotypes(N[k], poset, sampling_param=lambda_s, 
                                   lambdas=lambdas, eps=eps)
  
  Ne = foreach (i = 1:N[k], .combine='c') %dopar% {
    w = importance_weight(simulated_obs$obs_events[i, ], L=L, poset, lambdas,
                          lambda_s, eps, perturb_prob=perturb_prob,
                          sampling_time=NULL, sampling=sampling, seed=10)
    return(sum(w$w)^2 / sum(w$w^2))
  }
  effective_sample_size[[k]] = Ne
}
runtime = as.numeric(difftime(Sys.time(), t0, units='mins'))

###############################################################################
### TEST 6
###############################################################################
# Effective sample size as function of (i) compatible/incompatible status, 
# (ii) number of mutations (1's) per observation, (iii) number of mutations
# that could be added or (iv) removed, and (iv) hamming distance to compatible
# genotype that would be generated by either adding or removing

sampling = 'add-remove'
poset_trans_closed = trans_closure(poset)
parents = get_parents(poset_trans_closed)
childreen = get_childreen(poset_trans_closed)

set.seed(42)
n = 100
res = matrix(NA, nrow=n, ncol=9)
compatible = character(length=n)
for (i in 1:n) {
  genotype = simulated_obs$obs_events[i, ]
  res[i, 1] = sum(genotype == 1)
  compatible[i] = ifelse(is_compatible(genotype, poset), "C", "I")
  t0 = Sys.time()
  w = importance_weight(genotype=genotype, L=L, poset=poset, lambdas=lambdas,
                        lambda_s=lambda_s, eps=eps, sampling_time=NULL, 
                        perturb_prob=perturb_prob, sampling=sampling, seed=10)
  Ne = sum(w$w)^2 / sum(w$w^2)
  res[i, 2] = as.numeric(difftime(Sys.time(), t0, units="secs"))
  res[i, 3] = Ne*100/L
  
  # contrasting with known sampling time
  t0 = Sys.time()
  w = importance_weight(genotype=genotype, L=L, poset=poset,
                        lambdas=lambdas, lambda_s=lambda_s, eps=eps,
                        sampling_time=simulated_obs$T_sampling[i], 
                        perturb_prob=perturb_prob, sampling=sampling, seed=10)
  Ne = sum(w$w)^2 / sum(w$w^2)
  res[i, 4] = as.numeric(difftime(Sys.time(), t0, units="secs"))
  res[i, 5] = Ne*100/L
  
  compatible_by_adding = add_relations(genotype, childreen)
  compatible_by_removing = remove_relations(genotype, parents)
  dist_add = hamming_dist(genotype, compatible_by_adding)
  dist_remove = hamming_dist(genotype, compatible_by_removing)
  # Number of mutations on the compatible genotype potentially generated by 
  # adding mutations, and hamming distance
  res[i, 6] = sum(compatible_by_adding == 1)
  res[i, 7] = dist_add
  # Number of mutations on the compatible genotype potentially generated by 
  # removing mutations, and hamming distance
  res[i, 8] = sum(compatible_by_removing == 1)
  res[i, 9] = dist_remove
}

df = 
  data.frame("genotype"=
               apply(simulated_obs$obs_events[1:n, ], 1, paste, collapse=" "),
             "status"=compatible, "observed"=res[, 1], "runtime"=res[, 2],
             "Ne"=res[, 3], "runtime.Ts"=res[, 4], "Ne.Ts"=res[, 5],
             "observed.add"=res[, 6], "dist.add"=res[, 7],
             "observed.remove"=res[, 8], "dist.remove"=res[, 9])

###############################################################################
### TEST 7
###############################################################################
# MC sequential vs parallel execution 
library(doMC)
thrds = 4
registerDoMC(thrds)

set.seed(10)
t0 <- Sys.time()
res = MCEM_hcbn(poset, simulated_obs$obs_events, sampling_times=NULL, 
                lambda_s=lambda_s, max_iter=100, burn_in=0.8, L=1000, 
                sampling=sampling, parallel=FALSE, max_lambda_val=1e6, 
                verbose=TRUE)
runtime = as.numeric(difftime(Sys.time(), t0, units='mins'))
# [1] 0.1748256, p=5, N = 100, L = 100, iter = 100 (sequential V1), llhood=-959.3795, avg_llhood=-936.5417, eps=0.05618664, avg_eps=0.05370221
# [1] 0.1796839, p=5, N = 100, L = 100, iter = 100 (sequential V2), llhood=-959.3795, avg_llhood=-936.5417, eps=0.05618664, avg_eps=0.05370221
# [1] 0.3217676, p=5, N = 100, L = 100, iter = 100 (parallel), llhood=-956.5396, avg_llhood=-946.7611, eps=0.05619123, avg_eps=0.05347745
# [1] 1.702708, p=5, N = 1000, L = 100, iter = 100 (sequential V2), llhood=-9588.573, avg_llhood=-9672.281, eps=0.06164169, avg_eps=0.0621152
# [1] 1.898503, p=5, N = 1000, L = 100, iter = 100 (parallel), llhood=-9739.016, avg_llhood=-9716.249, eps=0.06312477, avg_eps=0.06283835
# [1] 0.266757, p=10, N = 100, L = 100, iter = 100 (sequential V2)
# [1] 0.3917821, p=10, N = 100, L = 100, iter = 100 (parallel), llhood=-3261.107, avg_llhood=-3169.789, eps=0.08171619, avg_eps=0.07664412
# [1] 0.9910405, p=5, N = 100, L = 1000, iter = 100 (sequential V2) llhood=-810.537, avg_llhood=-824.4606, eps=0.03436807, avg_eps=0.0358638
# [1] 0.7937994, p=5, N = 100, L = 1000, iter = 100 (parallel), llhood=-840.768, avg_llhood=-825.9859, eps=0.03763364, avg_eps=0.03598054
# [1] 9.935725, p=5, N = 1000, L = 1000, iter = 100 (sequential V2), llhood=-8540.092, avg_llhood=-8527.665, eps=0.04532772, avg_eps=0.04532164
# [1] 5.991118, p=5, N = 1000, L = 1000, iter = 100 (parallel), llhood=-8489.843, avg_llhood=-8514.486, eps=0.04481564, avg_eps=0.04506162

mean(abs(res$lambdas - lambdas))/mean(lambdas)
# [1] 0.06016701, p=5, N = 100, L = 100, iter = 100 (sequential V1)
# [1] 0.06016701, p=5, N = 100, L = 100, iter = 100 (sequential V2)
# [1] 0.05343616, p=5, N = 100, L = 100, iter = 100 (parallel)
# [1] 0.1243405, p=5, N = 1000, L = 100, iter = 100 (sequential V2)
# [1] 0.1122377, p=5, N = 1000, L = 100, iter = 100 (parallel)
# [1] 0.7943565, p=10, N = 100, L = 100, iter = 100 (sequential V2)
# [1] 0.736323, p=10, N = 100, L = 100, iter = 100 (parallel)
# [1] 0.06936015, p=5, N = 100, L = 1000, iter = 100 (sequential V2)
# [1] 0.0781487, p=5, N = 100, L = 1000, iter = 100 (parallel)
# [1] 0.07730737, p=5, N = 1000, L = 1000, iter = 100 (sequential V2)
# [1] 0.07993592, p=5, N = 1000, L = 1000, iter = 100 (parallel)

mean(abs(res$avg_lambdas - lambdas))/mean(lambdas)
# [1] 0.05160423, p=5, N = 100, L = 100, iter = 100 (sequential V1)
# [1] 0.05160423, p=5, N = 100, L = 100, iter = 100 (sequential V2)
# [1] 0.05316635, p=5, N = 100, L = 100, iter = 100 (parallel)
# [1] 0.1116793, p=5, N = 1000, L = 100, iter = 100 (sequential V2)
# [1] 0.1153271, p=5, N = 1000, L = 100, iter = 100 (parallel)
# [1] 0.6395589, p=10, N = 100, L = 100, iter = 100 (sequential V2)
# [1] 0.6401699, p=10, N = 100, L = 100, iter = 100 (parallel)
# [1] 0.07452685, p=5, N = 100, L = 1000, iter = 100 (sequential V2)
# [1] 0.07595211, p=5, N = 100, L = 1000, iter = 100 (parallel)
# [1] 0.07941715, p=5, N = 1000, L = 1000, iter = 100 (sequential V2)
# [1] 0.07881623, p=5, N = 1000, L = 1000, iter = 100 (parallel)