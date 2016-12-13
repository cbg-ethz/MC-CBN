library(mccbn)
library(doMC)
registerDoMC(4)

################## INPUT OPTIONS ##################
L = 10                              # number of repetitions
p = seq(4, 12, 2)                   # number of events
lambda_s = 1                        # sampling rate
N = unlist(lapply(50*p, min, 1000))  # number of observations / genotypes
eps = 0.05
##################################################

# Set seed for reproducibility
set.seed(47)

obs_log_likelihood = matrix(0, nrow=p, ncol=L)
relative_abs_error = matrix(0, nrow=p, ncol=L)
for (i in 1:length(p)) {
  res_list = foreach(j = 1:L, .combine=rbind) %dopar% {
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
}

