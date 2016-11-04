rm(list=ls())
library(graph)
library(mccbn)

set.seed(10)


##  we set the number of mutations and genotypes
# number of mutations
p = 10
# number of observed genotypes
N = 100

# generate a random poset with p nodes
true_poset = random_poset(p)

# plot the poset
plot_poset(true_poset)

## True mutation rates are randomly generated using a uniform distribution.
# Sequencing rate is set to one
lambda_s = 1

# generate p random mutation rates uniformly distributed between lambda_s/3 to 3lambda_s.  
lambdas = runif(p, 1/3*lambda_s, 3*lambda_s)


## N genotypes are simulated from the true CBN model i.e., poset and mutation rates using sample_genotypes
# Simulate genotypes and sequencing times consistent with poset and mutation rates
simGenotypes = sample_genotypes(N, true_poset, sampling_param = lambda_s, lambdas=lambdas)


# estimate mutation rates for a fixed poset 
est_lambda = estimate_mutation_rates(true_poset, simGenotypes$obs_events, simGenotypes$T_sampling) 

## We compare relative absoulte errors of estimates using the following command
abs(est_lambda$par - lambdas)/lambdas


# network learning on genotypes and sequencing times
# The following command takes around one minutes on a personal laptop
t1 = proc.time()
fit = learn_network(simGenotypes$obs_events, simGenotypes$T_sampling) 
proc.time()



## In order to obtain MLE poset(network), we postprocess the output of learn_network as follows:
  # MLE network
mle_index = which.max(fit$logliks)
plot_poset(fit$posets[[mle_index]])
