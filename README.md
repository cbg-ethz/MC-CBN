#MC-CBN
MC-CBN performs an efficient large-scale inference of the continuous time conjunctive Bayesian network (CBN). MC-CBN uses Monte Carlo expectation-maximization algorithm for inference of the CBN model.

##Usage
In this section, we explain how to use this R package on a simple example. In this example, we first generate a random poset as the true poset and assign to each node a mutation rate. Then, _N_ genotypes and corresponding sequencing times are simulated from this poset. Then, in the first task, we estimate the mutation rates for the given _true poset_. In the second task, we estimate the maximum likelihood (ML) poset for the observed data.

First, we load the necessary libraries and set the random seed for the reproducibility.   
```
library(graph)
library(mccbn)

set.seed(10)
```
Then we set the number of mutations and  genotypes
```
# number of mutations
p = 10
# number of observed genotypes
N = 100
```
We generate a random poset (partially ordered set) with _p_ nodes using the function _random_poset_ 
```
# generate a random poset with p nodes
true_poset = random_poset(p)
```
This poset is our _true poset_ in this example. Function _plot_poset_ can be used to plot this poset. 
```
# plot the poset
plot_poset(true_poset)
```
True mutation rates are randomly generated using a uniform distribution. 
```
# Sequencing rate is set to one
lambda_s = 1

# generate p random mutation rates uniformly distributed between lambda_s/3 to 3lambda_s.  
lambdas = runif(p, 1/3*lambda_s, 3*lambda_s)
```
_N_ genotypes are simulated from the true CBN model i.e., _poset_ and _mutation rates_ using _sample_genotypes_
```
# Simulate genotypes and sequencing times consistent with poset and mutation rates
simGenotypes = sample_genotypes(N, true_poset, sampling_param = lambda_s, lambdas=lambdas)
```
Now weestimate mutation rates using _estimate_mutation_rates_
```
# estimate mutation rates for a fixed poset 
est_lambda = estimate_mutation_rates(true_poset, simGenotypes$obs_events, simGenotypes$T_sampling) 
```
We simply compare relative absoulte errors of estimates using the following command
```
abs(est_lambda$par - lambdas)/lambdas
```
_learn_network_ can be used to perform poset learning using genotypes and sequencing times
```
# network learning on genotypes and sequencing times
# The following command takes around one minutes on a personal laptop
t1 = proc.time()
fit = learn_network(simGenotypes$obs_events, simGenotypes$T_sampling) 
proc.time()
```
In order to obtain MLE poset(network), we postprocess the output of _learn_network_ as follows:
```
# MLE network
mle_index = which.max(fit$logliks)
plot_poset(fit$posets[[mle_index]])
```

### Contributions
- [Hesam Montazeri](https://www.bsse.ethz.ch/cbg/group/people/person-detail.html?persid=168604)
- [Niko Beerenwinkel](http://www.bsse.ethz.ch/cbg/group/people/person-detail.html?persid=149417)

 ###Contact
```
Hesam Montazeri
hesam.montazeri (at) bsse.ethz.ch
```
