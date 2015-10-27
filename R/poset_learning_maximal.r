# 
# #' violation_freqs
# #' @export
# violation_freqs <- function(obs_events) {
#   p = ncol(obs_events)
#   N = nrow(obs_events)
#   violations <- matrix(0, p, p)
#   for(i in 1:p) {
#     for(j in 1:p) {
#       violations[i, j] = sum(obs_events[, i] == 0 & obs_events[, j] == 1)
#     }
#   }
#   diag(violations) = N+1
#   violations/N
# }


#' maximal_poset
#' @export
maximal_poset <- function(violations, eps) {
  p = ncol(violations)
  poset = matrix(0, p, p)
  for(i in 1:p) {
    for(j in 1:p) {
      if(violations[i, j] <= eps & poset[j, i] == 0) {
        
        parents_i = which(poset[, i] == 1)
        
        # check if adding (i,j) introduces a cycle                          
        if(all(poset[j, parents_i] == 0) ) 
        {
          poset[parents_i, j] = 1
          poset[i, j] = 1  
          
        } 
      }
    }
  }
  
  diag(poset) = 0
  trans_reduction(poset)
}


is_identical <- function(poset, posets) {
  for( tmp_poset in posets) {
    if(identical(tmp_poset, poset) ) {
      return(TRUE)
    }
  }
  return(FALSE)
}

# 
# candidate_posets <- function( obs_events ) {
#   posets <- list()
#   violations = violation_freqs(obs_events)
#   epsilons = sort(unique(c(0, violations)))
#   print(length(epsilons))  
#   for(eps in epsilons) {
#     poset = maximal_poset(violations, eps)
#     if( is_identical(poset, posets) == FALSE) {
#       posets[[length(posets) + 1]] = poset
#     }
#   }
#   
#   print(length(posets))
#   posets
# }


violation_freqs_w <- function (obs_events, obs_weights) {
  p = ncol(obs_events)
  N = sum(obs_weights)
  violations <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      violations[i, j] = sum(obs_weights[obs_events[, i] == 0 & obs_events[, j] == 1])
    }
  }
  diag(violations) = N + 1
  violations/N
}



candidate_posets <- function (obs_events, obs_weights) {
  posets <- list()
  violations = violation_freqs_w(obs_events, obs_weights)
  epsilons = sort(unique(c(0, violations)))
  print(length(epsilons))
  for (eps in epsilons) {
    poset = maximal_poset(violations, eps)
    if (is_identical(poset, posets) == FALSE) {
      posets[[length(posets) + 1]] = poset
    }
  }
  print(length(posets))
  posets
}

#' discretize_times
#' @export 
discretize_times <- function(times_, D) {
  N = length(times_)
  if(is.na(D)) {
    D = ( max(times_) - min(times_) ) /50  
  }
  
  t_indexes = rep(0, N)
  
  times2 = D * round(times_/D)
  quan_times = unique(times2)
  for(i in 1:N) {
    t_indexes[which(times2 == quan_times[i])] = i  
  }
  
  list(times = quan_times, t_indexes=t_indexes, all_times = quan_times[t_indexes])
}

all_maximal_posets <- function(obs_events, sampling_times, D=NA, par=TRUE, num_par_posets=10, optim_method="minqa", control=list())
{
  if(ncol(obs_events) > 19) {
    stop("Error, the exact method is only applicable to posets with less than 20 mutations!")
  }
  posets = candidate_posets( obs_events )
  
  nr_posets = length(posets)
  
  # initialize the variables
  alphas= logliks = rep(0, nr_posets)
  lambdas_mat = matrix(0, nr_posets, ncol(obs_events))
  fits = list()
  
  if(par) {
    mcoptions <- list(cores=num_par_posets)  
    results = foreach(i = 1:nr_posets, .options.multicore = mcoptions) %dopar%
    {
      poset = posets[[i]]
      compatible_geno = compatible_genotypes(obs_events, poset)
      print(poset)
      t1 = proc.time()
      
      dtimes = discretize_times(sampling_times[compatible_geno$compatible_indexes], D)
      fit = estimate_lambda(obs_events[compatible_geno$compatible_indexes,], 
                            dtimes, poset,
                            init_lambda = NA, grad = NULL, 
                            likelihood_fn=likelihood, maxfun=maxfun, optim_method=optim_method)
      
      lambdas = fit$par
      print(fit)
      print(proc.time() - t1)
      t1 = proc.time()
#       cur_loglik = loglikelihood_cbn(poset, lambdas, sampling_times, obs_events, D)
      ll_control = list(ll_method = "exact", D=D)
      cur_loglik = loglike_mixture_model(poset, lambdas, obs_events, sampling_times, ll_control) 
      print(cur_loglik)
      print(proc.time() - t1)
      print(paste("*** done ", sep = ""))
      
      list(fit=fit, alpha=cur_loglik$alpha, ll=cur_loglik$ll, lambdas=lambdas)
    }
      # update the variables for the poset i
    for(i in 1:nr_posets) {
      fits[[i]] = results[[i]]$fit
      alphas[i] = results[[i]]$alpha
      logliks[i] = results[[i]]$ll
      lambdas_mat[i, ] = results[[i]]$lambdas
    }
  } else {
    print("Serial computation: argument num_par_posets is ignored!")
    for(i in 1:nr_posets) # for debugging simple for has to be used
    {
      poset = posets[[i]]
      compatible_geno = compatible_genotypes(obs_events, poset)
      print(poset)
      t1 = proc.time()
      dtimes = discretize_times(sampling_times[compatible_geno$compatible_indexes], D)
      
#       fit = fit_params(poset, obs_events[compatible_geno$compatible_indexes,], times, err_model="no-error", D=D, optim_method = optim_method)
# 
# estimate_lambda <- function(obs_events, dtimes, poset_, likelihood_fn = likelihood, init_lambda = NA, init_eps=0.3, grad=NULL, 
#                             optim_method="minqa", control=list(), with_eps=FALSE, ...) 
  
  fit = estimate_lambda(obs_events[compatible_geno$compatible_indexes,], 
                        dtimes, poset, likelihood_fn=likelihood, 
                        init_lambda = NA, grad = NULL, 
                       optim_method=optim_method,  control=control)


#        fit = estimate_lambda(obs_events[compatible_geno$compatible_indexes,], 
#                             dtimes, poset,
#                             init_lambda = NA, grad = NULL, 
#                             likelihood_fn=likelihood, maxfun=maxfun, optim_method=optim_method)
      lambdas = fit$par
      print(fit)
      print(proc.time() - t1)
      t1 = proc.time()
#       loglike_mixture_model <- function(poset, lambdas, obs_events, times, D)
#       cur_loglik = loglikelihood_cbn(poset, lambdas, sampling_times, obs_events, D)
      ll_control = list(ll_method = "exact", D=D)
      cur_loglik = loglike_mixture_model(poset, lambdas, obs_events, sampling_times, ll_control) 
      print(cur_loglik)
      print(proc.time() - t1)
      print(paste("*** done ", sep = ""))
      
      # update the variables for the poset i
      fits[[i]] = fit
      alphas[i] = cur_loglik$alpha
      logliks[i] = cur_loglik$ll
      lambdas_mat[i, ] = lambdas
    }
  }

list(posets = posets, alphas = alphas, fits = fits,
     logliks = logliks, lambdas_mat = lambdas_mat, obs_events = obs_events,
     sampling_times = sampling_times)
}
