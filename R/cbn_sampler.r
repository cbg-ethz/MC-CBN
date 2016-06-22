
#' make_linear_poset
#' @param p the number of events
#' @export
make_linear_poset <- function(p) {
  Poset <- matrix(0, p, p)
  if(p > 1){
    for(i in 1:(p-1))
      Poset[i, i+1] = 1
  }
  Poset
}

#' Making an empty poset
#' 
#' @description This function creates an empty poset with p nodes and no edges.
#' @param p the number of events
#' @return the empty poset of size \code{p}
#' @examples 
#' make_empty_poset(4)
#' @export 
make_empty_poset <- function(p) {
  matrix(0, p, p)
}


#' make_random_poset
#' @param p the number of events
#' @export
make_random_poset <- function(p) {
  poset <- matrix(0, p, p)
  if(p == 4){
    poset[1, 3] = poset[2, 3] = poset[2, 4] = 1
  }
  
  if(p == 6) {
    poset[1, 4] = poset[2, 4] = poset[3, 4] = 1
    poset[3, 5] = 1
    poset[4, 6] = poset[5, 6] = 1
  }
  poset
}


#' sampling.expo
#' @export
sampling.expo <- function(n, mean_s, mut_times=NULL) {
  rexp(n, 1/mean_s)
}

#' sampling.norm
#' @export
sampling.norm <- function(n, mean_s, mut_times=NULL) {
  #rnorm(n, mean_s, 0.1*(abs(mean_s)+1) )
  if(mean_s <= 0) {
    print("Error, mean_s should be positive!")
  }
  rtruncnorm(n, a=0, b=Inf, mean=mean_s, sd=0.1*mean_s )
}

#' sampling.const 
#' @export
sampling.const <- function(n, mean_s, mut_times=NULL) {
  rep(mean_s, n)
}

#' sampling.const 
#' @export
sampling.unif <- function(n, lambdas, mut_times=NULL) {
  max = sum(lambdas/lambdas)
  runif(n, max/100, max )
}


#' sampling.dep 
#' @export
sampling.dep <- function(n, random_fraction, mut_times) {
  if(is.null(mut_times) )
    return(NULL)
  
  returned_sampling_times = t(apply(mut_times, 1, function(x) { c(max(x[1:2]), runif(1, min(x)/1.1  ,max(x)*1.1) )} ) )
  
  component_selector_idx = sample.int(2, n, replace=TRUE, prob=c(1-random_fraction, random_fraction) ) 
  
  final_sampling_T = rep(0, n)
  final_sampling_T[component_selector_idx==1] = returned_sampling_times[component_selector_idx==1, 1]
  final_sampling_T[component_selector_idx==2] = returned_sampling_times[component_selector_idx==2, 2]
  final_sampling_T
}

#' sample_timed_genotypes_with_eps
#' @export
sample_genotypes <- function (n, poset, sampling_param, lambdas, sampling_fn=sampling.expo, eps=0.0) 
{
  p = length(lambdas)
  T_events <- matrix(0, n, p)
  T_sampling <- sampling_fn(n, sampling_param, NULL)
  for (i in 1:p) {
    T_events[, i] <- rexp(n, lambdas[i])
  }
  T_sum_events <- matrix(0, n, p)
  topo_path = my.topological.sort(poset)
  for (e in topo_path) {
    parents <- which(poset[, e] == 1)
    if (length(parents) == 0) {
      T_sum_events[, e] = T_events[, e]
    }
    else if (length(parents) == 1) {
      T_sum_events[, e] = T_events[, e] + T_sum_events[, 
                                                       parents]
    }
    else {
      T_sum_events[, e] = T_events[, e] + apply(T_sum_events[, 
                                                             parents], 1, max)
    }
  }
  if (is.null(T_sampling)) 
    T_sampling <- sampling_fn(n, sampling_param, T_sum_events)
  obs_events <- matrix(0, n, p)
  for (i in 1:p) {
    obs_events[, i] <- rbinom(n, as.numeric(T_sum_events[, i] <= T_sampling), 1-eps)
  }
  genotype_list = apply(obs_events, 1, function(x) {
    which(x == 1)
  })
  list(n = n, p = p, T_sampling = T_sampling, obs_events = obs_events, 
       genotype_list = genotype_list, sampling_param = sampling_param, 
       lambdas = lambdas, T_events = T_events, T_sum_events = T_sum_events, eps=eps)
}