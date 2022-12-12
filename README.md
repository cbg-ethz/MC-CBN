# MC-CBN
MC-CBN performs an efficient large-scale inference of the continuous time conjunctive Bayesian network (CBN). MC-CBN uses Monte Carlo expectation-maximization algorithm for inference of the CBN model.

## Usage
In this section, we explain how to use the R package on a simple example. In this example, we first generate a random poset as the true poset and assign to each node a mutation rate. Then, _N_ genotypes and corresponding sequencing times are simulated from this poset. Then, in the first task, we estimate the mutation rates for the given _true poset_. In the second task, we estimate the maximum likelihood (ML) poset for the observed data.

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

### OT-CBN
_N_ genotypes are simulated from the true CBN model i.e., _poset_ and _mutation rates_ using _sample_genotypes_
```
# Simulate genotypes and sequencing times consistent with poset and mutation rates
simGenotypes = sample_genotypes(N, true_poset, sampling_param = lambda_s, lambdas=lambdas)
```
Now we estimate mutation rates using _estimate_mutation_rates_
```
# estimate mutation rates for a fixed poset 
est_lambda = estimate_mutation_rates(true_poset, simGenotypes$obs_events, simGenotypes$T_sampling) 
```
We simply compare relative absoulte errors of estimates using the following command
```
abs(est_lambda$lambda - lambdas)/lambdas
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
### H-CBN2
Now, we generate genotypes from the true CBN model, but in addition, we simulate noisy observations
```
simGenotypes = sample_genotypes(N, true_poset, sampling_param = lambda_s, lambdas=lambdas, eps=0.01)
```

We estimate model parameters using _MCEM.hcbn_
```
lambda0 = runif(p, 1/3*lambda_s, 3*lambda_s)
est_params = MCEM.hcbn(lambda0, true_poset, simGenotypes$obs_events, lambda.s=lambda_s, L=100, eps=0.01)
```

To learn the poset from the data, we use _adaptive.simulated.annealing_
```
# Initial poset
posets = candidate_posets(simGenotypes$obs_events, rep(1, N), 0.9)
poset0 = posets[[length(posets)]]

fit = adaptive.simulated.annealing(poset0, simGenotypes$obs_events, L=100,
				   max.iter.asa=10, seed=10L)
```

## Installation

### Dependencies

MC-CBN requires the following software:

- C++ compiler with support for c++11 and OpenMP
- Boost C++ library (tested with version 1.74.0)
- (optional) Intel Math Kernel Library (MKL)

and R-packages for Bioconductor:
- graph
- Rgraphviz

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version="3.16")
BiocManager::install(c("graph", "Rgraphviz"))
```

To install MC-CBN, you can use the latest tarball release, as well as R-package `remotes` to install dependencies, e.g:
```
install.packages("remotes")
remotes::install_url("https://github.com/cbg-ethz/MC-CBN/releases/download/v2.1.0/mccbn_2.1.0.tar.gz", dependencies=TRUE)
```


To compile with MKL, necessary compiler and linker flags can be determined with the help of the Intel MKL [Link Line Advisor](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl/link-line-advisor.html), e.g.:

```
remotes::install_url("https://github.com/cbg-ethz/MC-CBN/releases/download/v2.1.0/mccbn_2.1.0.tar.gz", dependencies=TRUE,
                     configure.args="--with-mklcxxflags=\"-DMKL_ILP64 -m64 -I${MKLROOT}/include\" --with-mklldflags=\"-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl\"")
```

### Installation from source

Requirements:
- Autotools (Autoconf >= 2.69, Automake >= 1.15, Automake-archive)
- R-package `remotes`

After cloning the repository, go to the corresponding directory and type:

```
autoreconf -vif
```

Then build the package:
```
R CMD build .
````

Finally, install the package and other dependencies using R-package `remotes`:
```
install.packages("remotes")
remotes::install_local("mccbn_<version>.tar.gz", dependencies=TRUE)
```

Once again to compile with MKL, the appropriate flags can be passed using `--configure-args`, e.g:

```
remotes::install_local("mccbn_<version>.tar.gz", dependencies=TRUE,
                       configure.args="--with-mklcxxflags=\"-DMKL_ILP64 -m64 -I${MKLROOT}/include\" --with-mklldflags=\"-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl\"")
```

### Contributions
- [Hesam Montazeri](https://www.bsse.ethz.ch/cbg/group/people/person-detail.html?persid=168604)
- [Susana Posada Cespedes](https://www.bsse.ethz.ch/cbg/group/people/person-detail.html?persid=192769)
- [Niko Beerenwinkel](http://www.bsse.ethz.ch/cbg/group/people/person-detail.html?persid=149417)

### Contact
```
Susana Posada Cespedes
susana.posada (at) bsse.ethz.ch
```
