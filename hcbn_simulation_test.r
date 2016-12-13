library(mccbn)
library(doMC)
registerDoMC(4)

################## INPUT OPTIONS ##################
L = 3                               # number of repetitions
p = c(4, 6) #seq(4, 12, 2)                   # number of events
lambda_s = 1                        # sampling rate
N = unlist(lapply(50*p, min, 1000))  # number of observations / genotypes
eps = 0.05
mccbn_path = "/Users/susanap/Documents/software/MC-CBN"
##################################################

source(file.path(mccbn_path, "hcbn_functions.r"))

# Set seed for reproducibility
set.seed(47)

obs_llhood = matrix(0, nrow=length(p), ncol=L)
relative_abs_error = matrix(0, nrow=length(p), ncol=L)
for (i in 1:length(p)) {
  res = foreach(j = 1:L, .combine=rbind) %dopar% {
    poset = make_random_poset(p[i])
    lambdas = runif(p[i], 1/3*lambda_s, 3*lambda_s)
    simulated_obs = sample_genotypes(N[i], poset, sampling_param=lambda_s, 
                                     lambdas=lambdas, eps=eps)
    ret = MCMC_hcbn(poset, simulated_obs$obs_events)
    llhood = obs_log_likelihood(simulated_obs$obs_events, poset, lambdas,
                                lambda_s, eps, exact=TRUE)
    error = mean(abs(ret$lambdas - lambdas))/mean(lambdas)
    return(c(llhood, error))
  }
  obs_llhood[i, ] = res[, 1]
  relative_abs_error[i, ] = res[, 2]
}

