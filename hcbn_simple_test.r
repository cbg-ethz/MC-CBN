library(mccbn)

# Set seed for reproducibility
set.seed(42)

# input
p = 20                     # number of events
poset = random_poset(p)    # true poset
lambda_s = 1               # sampling rate
N = 1000                   # number of observations / genotypes
L = 1000                   # number of samples
eps = 0.05
mccbn_path = "/Users/susanap/Documents/software/MC-CBN"

# Load functions 
source(file.path(mccbn_path, "hcbn_functions.r"))

# simulate data set 
lambdas = runif(p, 1/3*lambda_s, 3*lambda_s)
simulated_obs = sample_genotypes(N, poset, sampling_param=lambda_s, lambdas=lambdas,
                                 eps=eps)
# Effective sample size
idx = 2
t0 = Sys.time()
w = importance_weight(genotype=simulated_obs$obs_events[idx, ], L=L, poset=poset,
                      lambdas=lambdas, lambda_s=lambda_s, eps=eps,
                      sampling_time=NULL, sampling='add-remove')
Ne = sum(w$w)^2 / sum(w$w^2)
print(difftime(Sys.time(), t0, units='mins'))
print(Ne)
