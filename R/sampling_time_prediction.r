.tlikelihood <- function(t, lambdas, m) {
  P = pexp(t, lambdas)
  
  logP = m*P + (1-m)*(1-P)
  logP[ logP==0 ] = 10^-50
  -sum( log( logP ) )
}

.fit_exp <- function(x, lambda_s=1, verbose=TRUE) {
  M = sum(x == 1)
  N = length(x)
  
  if(M == 0 || N==M) {
    stop("Error! M == 0 or N==M. Get out of my sight!")
  }
  
  lambda = M * lambda_s / (N-M)
  expected = N* (lambda/(lambda+lambda_s) )
  list(lambda=lambda,expected = expected)
}


.MLE_sampling_time_for_genotype <- function(x, rates, itime=1) {
  optim(itime, .tlikelihood,  method = "L-BFGS-B", lower = 0.000001, upper=50, lambdas=rates, m=x)$par
}

estimate_sampling_times <- function(X, N) {
  # we assume a naive Bayes model. We first esimate the rate of each mutation.
  rates = apply(X, 2, function(x) { .fit_exp(x, lambda_s=1, verbose=TRUE)$lambda}) 
  
  # Then we find the ML estimation of sampling time for each genotype with the exponential rates
  ranges = t(apply(X, 1, expected_sampling_time_interval_for_genotype, rates=rates, N))
  
  apply(ranges, 1, function(x) {
    if(is.infinite(x[2]) ) {
      return(1.1*x[1] )
    } else{
      return(mean(x))
    }
  })
}

expected_sampling_time_interval_for_genotype <- function(geno_, rates, N=1000) {
  geno = which(geno_ == 1)
  T = c()
  for(rate in rates) {
    T = cbind(T, rexp(N, rate) )
  }
  
  if(length(geno) == 0) {
    T_obs = 0
    T_unobs = mean(apply(T, 1, min)  )  
  } else {
    T_obs = mean(apply(T, 1, function(x){ max(x[geno]) }) )  
    
    if(length(geno) != length(rates)) {
      T_unobs = mean(apply(T, 1, function(x){ min(x[-geno]) }) )  
    } else  {
      T_unobs = Inf
    }
      
  }
#     
#   print(T_obs)
#   print(T_unobs)
#   print(T_obs + T_unobs)
#   
#   print( )
  c(T_obs, T_obs + T_unobs)
}


estimate_sampling_times2 <- function(X) {
  # we assume a naive Bayes model. We first esimate the rate of each mutation.
  rates = apply(X, 2, function(x) { .fit_exp(x, lambda_s=1, verbose=TRUE)$lambda}) 
  
  # Then we find the ML estimation of sampling time for each genotype with the exponential rates
  apply(X, 1, .MLE_sampling_time_for_genotype, rates=rates, itime=1)
}



################### Expected sampling time by importance sampling ######################

expected_sampling_time_importance_sampling <- function(genotype, poset, lambda, nrOfSamples ) {
  topo_path = my.topological.sort(poset) -1
  stimes = rexp(nrOfSamples, 1)
  
  T = list(T=matrix(0, 0, ncol(poset)), densities=c())
  for(stime in stimes) {
    t_ = .Call("drawHiddenVarsSamples", 1, genotype, stime, lambda,  poset, topo_path)  
    
    T$T  = rbind(T$T, t_$T)
    T$densities = c(T$densities, t_$densities)
  }
  
  values = cbn_density_(T$T, lambda) - T$densities
  #   logp = logsumexp(values)$logR
  
  #    weights = exp(cbn_density_(T$T, lambda) - T$densities)
  
  #   (sum(stimes * weights)/ nrOfSamples) /exp(logp) 
  #     sum(stimes * weights)/ sum(weights)
  
  logp = logsumexp(values)$logR
  exp(logsumexp(log(stimes) + values)$logR - logp)
}
