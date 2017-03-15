library(mccbn)

# Set seed for reproducibility
set.seed(42)

# input
p = 8                      # number of events
poset = random_poset(p)    # true poset
lambda_s = 1               # sampling rate
N = 1000                   # number of observations / genotypes
L = 1000                   # number of samples
eps = 0.05
perturb_prob = 0
# Enter the path or set it as the working directory
mccbn_path = "/Users/susanap/Documents/software/MC-CBN"

# Load functions 
if (file.exists(file.path(mccbn_path, "hcbn_functions.r"))) {
  source(file.path(mccbn_path, "hcbn_functions.r"))
} else if (file.exists("hcbn_functions.r")) {
  source("hcbn_functions.r")
}


# simulate data set 
lambdas = runif(p, 1/3*lambda_s, 3*lambda_s)
simulated_obs = sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas,
                                 eps=eps)
###############################################################################
## Test 1 - Exploratory
## Effective sample size for particular observations and on average for 
## 100 obervations
###############################################################################
idx = 2
genotype = simulated_obs$obs_events[idx, ]
t0 = Sys.time()
w = importance_weight(genotype=genotype, L=L, poset=poset, lambdas=lambdas,
                      lambda_s=lambda_s, eps=eps, sampling_time=NULL, 
                      perturb_prob=perturb_prob, sampling='add-remove')
Ne = sum(w$w)^2 / sum(w$w^2)
print(difftime(Sys.time(), t0)) # Time difference of 2.484901 secs
print(Ne) 
# [1] 545.4277 (p=8, perturb_prob = 0, eps = 0.05)
# [1] 49.0073 (p=10, perturb_prob = 0.8, eps = 0.05) - WT
# [1] 115.4417 (p=10, perturb_prob = 0.8, eps = 0.05)
# [1] 150.1037 (p=10, perturb_prob = 0.8, eps = 0.1)

# contrasting with known sampling time
t0 = Sys.time()
w = importance_weight(genotype=genotype, L=L, poset=poset,
                      lambdas=lambdas, lambda_s=lambda_s, eps=eps,
                      sampling_time=simulated_obs$T_sampling[idx], 
                      perturb_prob=perturb_prob, sampling='add-remove')
Ne = sum(w$w)^2 / sum(w$w^2)
print(difftime(Sys.time(), t0)) # Time difference of 2.492921 secs
print(Ne) 
# [1] 974.5802 (p=8, perturb_prob = 0, eps = 0.05)
# [1] 238.5888 (p=10, perturb_prob = 0.8, eps = 0.05) - WT
# [1] 345.9204 (p=10, perturb_prob = 0.8, eps = 0.1)
# [1] 185.8936 (perturb_prob = 0.8, eps = 0.05) 

# Over more observations
t0 = Sys.time()
Ne = numeric(100)
for (i in 1:100) {
  w = importance_weight(genotype=simulated_obs$obs_events[i, ], L=L, poset=poset,
                        lambdas=lambdas, lambda_s=lambda_s, eps=eps, 
                        perturb_prob=perturb_prob, sampling_time=NULL, 
                        sampling='add-remove')
  Ne[i] = sum(w$w)^2 / sum(w$w^2)
}
cat(difftime(Sys.time(), t0, units='mins'), mean(Ne*100/L), sd(Ne*100/L), 
    max(Ne*100/L), min(Ne*100/L), sep=" ")

t0 = Sys.time()
Ne = numeric(100)
for (i in 1:100) {
  w = importance_weight(genotype=simulated_obs$obs_events[i, ], L=L, poset=poset,
                        lambdas=lambdas, lambda_s=lambda_s, eps=eps, 
                        perturb_prob=perturb_prob, 
                        sampling_time=simulated_obs$T_sampling[i], sampling='add-remove')
  Ne[i] = sum(w$w)^2 / sum(w$w^2)
}
cat(difftime(Sys.time(), t0, units='mins'), mean(Ne*100/L), sd(Ne*100/L), 
    max(Ne*100/L), min(Ne*100/L), sep=" ")
# [1] 244.9589 (p=10, perturb_prob = 0.8, eps = 0.05) - WT
# [1] 345.9204 (p=10, perturb_prob = 0.8, eps = 0.1)
# [1] 101.2962 (perturb_prob = 0.8, eps = 0.05) 

###############################################################################
## Test 2
## Effective sample size for many observations
###############################################################################

poset_trans_closed = trans_closure(poset)
parents = get_parents(poset_trans_closed)
childreen = get_childreen(poset_trans_closed)

set.seed(42)
n = 100
res = matrix(NA, nrow=n, ncol=9)
compatible = character(length=10)
for (i in 1:n) {
  genotype = simulated_obs$obs_events[i, ]
  res[i, 1] = sum(genotype == 1)
  compatible[i] = ifelse(is_compatible(genotype, poset), "C", "I")
  t0 = Sys.time()
  w = importance_weight(genotype=genotype, L=L, poset=poset, lambdas=lambdas,
                        lambda_s=lambda_s, eps=eps, sampling_time=NULL, 
                        perturb_prob=perturb_prob, sampling='add-remove')
  Ne = sum(w$w)^2 / sum(w$w^2)
  res[i, 2] = as.numeric(difftime(Sys.time(), t0, units="secs"))
  res[i, 3] = Ne*100/L
  
  # contrasting with known sampling time
  t0 = Sys.time()
  w = importance_weight(genotype=genotype, L=L, poset=poset,
                        lambdas=lambdas, lambda_s=lambda_s, eps=eps,
                        sampling_time=simulated_obs$T_sampling[i], 
                        perturb_prob=perturb_prob, sampling='add-remove')
  Ne = sum(w$w)^2 / sum(w$w^2)
  res[i, 4] = as.numeric(difftime(Sys.time(), t0, units="secs"))
  res[i, 5] = Ne*100/L
  
  compatible_by_adding = add_relations(genotype, childreen)
  compatible_by_removing = remove_relations(genotype, parents)
  dist_add = hamming_dist(genotype, compatible_by_adding)
  dist_remove = hamming_dist(genotype, compatible_by_removing)
  res[i, 6] = sum(compatible_by_adding == 1)
  res[i, 7] = dist_add
  res[i, 8] = sum(compatible_by_removing == 1)
  res[i, 9] = dist_remove
}

df = data.frame("genotype"=apply(simulated_obs$obs_events[1:n, ], 1, paste, collapse=" "),
                "status"=compatible, "observed"=res[, 1], "runtime"=res[, 2], 
                "Ne"=res[, 3], "runtime.Ts"=res[, 4], "Ne.Ts"=res[, 5], 
                "observed.add"=res[, 6], "dist.add"=res[, 7], 
                "observed.remove"=res[, 8], "dist.remove"=res[, 9])


