earn_network_boot <- function(obs_events, sampling_times, B = 50, weights=NULL, max_iter=200, zeta = 0.2, nrOfSamplesForEStep=50,
                               nrOfSamplesForLL = 100,   noise_model="empty", verbose=FALSE, min_compatible_geno_fraction = 0.5) {
  
  est_poset = est_poset_trans = matrix(0, ncol(obs_events), ncol(obs_events))
  N = nrow(obs_events)
  if(is.null(weights)) {
    weights = rep(1, nrow(obs_events))
  }
  
  for(i in 1:B) {
    
    indexes = sample(nrow(obs_events),  N, replace =TRUE )
    obs_events2 = obs_events[indexes, ]
    times2 = sampling_times[indexes]
    weights2 = weights[indexes]

    result = learn_network(obs_events2, times2, weights2, max_iter, zeta, nrOfSamplesForEStep,
                           nrOfSamplesForLL, noise_model, verbose, min_compatible_geno_fraction)   
    ml_index = which.max(result$loglik)
    mle_poset = result$posets[[ml_index]]
    
    est_poset = est_poset + mle_poset
    est_poset_trans = est_poset_trans + trans_closure(mle_poset)
  }
  list(poset= est_poset/B, poset_trans = est_poset_trans/B)
}

learn_network <- function(obs_events, sampling_times = NULL, weights=NULL, max_iter=100, zeta = 0.2, L=5,
                          nrOfSamplesForLL = 100,   noise_model="empty", verbose=FALSE, min_compatible_geno_fraction=0.5, maxLambdaValue=10^6, lambda_s=1.0) {
  
  sampling_times_available = is.null(sampling_times) == FALSE
  
  learn_network_internal(obs_events, sampling_times, max_iter, zeta, L, nrOfSamplesForLL, weights, removeZeroWeights=TRUE,   noise_model, 
      verbose, min_compatible_geno_fraction, maxLambdaValue=maxLambdaValue, lambda_s=lambda_s, sampling_times_available=sampling_times_available)  
}


learn_network_internal <- function(obs_events, sampling_times, max_iter=200, zeta = 0.2, nrOfSamplesForEStep=50,
                                    nrOfSamplesForLL = 100, weights=NULL, removeZeroWeights=TRUE, noise_model="empty", verbose=FALSE, 
                                    min_compatible_geno_fraction, maxLambdaValue, lambda_s, sampling_times_available)
{
  if(sampling_times_available == FALSE) {
    sampling_times = rep(0, nrow(obs_events))
  } 
  
  if(is.matrix(obs_events) == FALSE) {
    stop("Function all_maximal_posets_mcem: obs_events must be of 'matrix' type!")
  }
  
  if(is.null(weights)) {
    weights = rep(1, nrow(obs_events))
  }
  
  #### select observations with nonzero weights
  if(removeZeroWeights) {
    indexes = which(weights > 0)
    if(length(indexes) == 0){
      stop("No observation with nonzero weight!")
    } 
    
    obs_events = obs_events[indexes, , drop=FALSE]
    weights = weights[indexes]
    sampling_times = sampling_times[indexes]
  }
  
  posets = candidate_posets( obs_events, weights, min_compatible_geno_fraction )
  nr_posets = length(posets)
  
  # initialize the variables
  alphas= logliks = rep(0, nr_posets)
  lambdas_mat = matrix(0, nr_posets, ncol(obs_events))
  fits = list()
  
  max_loglike = -Inf # loglike is always negative
  
  for(i in 1:nr_posets) 
  {
    poset = posets[[i]]
    compatible_geno = compatible_genotypes(obs_events, poset)
    
    if(verbose) {
      print("**********************")
      print(poset)
    }
    
    t1 = proc.time()
    ilambda = NULL
    fit = NA
    lambdas = rep(NA, ncol(poset))
    
    if(noise_model == "empty"){
      nc_indexes = setdiff(1:nrow(obs_events), compatible_geno$compatible_indexes)
      
      geno_prob_noise = NA
      if(length(nc_indexes) > 0) {
        geno_prob_noise = genotype_probs_empty_poset(obs_events[nc_indexes, , drop=FALSE], sampling_times[nc_indexes], weights[nc_indexes], 
                                                     max_iter, zeta, nrOfSamplesForEStep, verbose, maxLambdaValue, lambda_s, sampling_times_available)
        
      }
    }
    
    logkile_incompatible = incompatible_loglike(poset, obs_events, sampling_times, weights, compatible_geno, noise_model, geno_prob_noise)
    
    if( (i > 1 && (max_loglike >  logkile_incompatible) ) || ( length(compatible_geno$compatible_indexes) <= 0 ) ) {
      cur_loglik = logkile_incompatible
      if(verbose==TRUE) {
        print("No parameter estimation for this poset. Upper bound likelihood of this poset is less than current maximum likelihood poset!")        
      }
    } else{
      fit = MCEM(poset, obs_events[compatible_geno$compatible_indexes, , drop=FALSE ], sampling_times[compatible_geno$compatible_indexes], max_iter=max_iter,
                 weights=weights[compatible_geno$compatible_indexes], zeta = zeta, ilambda=ilambda, nrOfSamples=nrOfSamplesForEStep, verbose=verbose, 
                 maxLambdaValue=maxLambdaValue, lambda_s=lambda_s, sampling_times_available=sampling_times_available)
      
      lambdas = fit$par
      if(verbose) {
        print(fit)
        print(proc.time() - t1)
      }
      t1 = proc.time()
      cur_loglik = loglike_mixture_model(poset, lambdas, obs_events, sampling_times, weights, nrOfSamplesForLL, compatible_geno, 
                                         logkile_incompatible, lambda_s, sampling_times_available) 
    }
    
    if(verbose) {
      print(cur_loglik)
      print(proc.time() - t1)
      print(paste("*** done ", sep = ""))
    }
    
    # update the variables for the poset i
    fits[[i]] = fit
    alphas[i] = cur_loglik$alpha
    logliks[i] = cur_loglik$ll
    lambdas_mat[i, ] = lambdas
    
    if(logliks[i]  > max_loglike) {
      max_loglike = logliks[i]
    }
  }
  
  
  list(posets = posets, alphas = alphas, fits = fits,
       logliks = logliks, lambdas_mat = lambdas_mat, obs_events = obs_events,
       sampling_times = sampling_times)
}



genotype_probs_empty_poset <- function (obs_events, sampling_times, weights, max_iter, zeta, nrOfSamplesForEStep, verbose,
                                        maxLambdaValue, lambda_s, sampling_times_available) {
  # all the observed genotypes.
  poset_noise = make_empty_poset(ncol(obs_events) ) 
  
  # 1) estimate the rates from all the data use MCEM function
  fit = MCEM(poset_noise, obs_events, sampling_times, max_iter=max_iter,
             weights=weights, zeta=zeta, ilambda=NULL, nrOfSamples=nrOfSamplesForEStep, verbose=verbose,
             maxLambdaValue=maxLambdaValue, lambda_s=lambda_s, sampling_times_available=sampling_times_available)
  lambdas_noise = fit$par
  
  geno_prob_noise = all_genotype_prob_for_empty_poset(lambdas_noise, obs_events, sampling_times, lambda_s, sampling_times_available,log.p=TRUE)
  geno_prob_noise 
}


genoLatticeSize_fast <- function(poset, already_transitive=FALSE) {
  if(already_transitive == FALSE) {
    poset = trans_closure(poset)  
  }
  .Call("genoLatticeSize", poset, FALSE) 
}


genoLatticeSize <- function(poset, already_transitive=FALSE) {
  if(already_transitive == FALSE) {
    poset = trans_closure(poset)  
  }
  genoLatticeSizeRecursion(poset)
}


genoLatticeSizeRecursion <- function(poset) {
  
  if(all(poset == 0)) {
    return(ncol(poset) * log(2) )
  }
  
  indegs = apply(poset, 2, sum)
  next_node = which.min( abs( indegs - mean(indegs) ) )[1]
  next_node = 1
  
  log_size1 = genoLatticeSizeRecursion(as.matrix(poset[-next_node,-next_node]) ) 
  # just in comment size1= exp(log_size1)
  
  indexes = unique(c(next_node, which(poset[next_node,] == 1), which(poset[,next_node] == 1)) )
  
  if(length(indexes) == 1) {
    # 2 * size1
    return(log(2) + log_size1)    
  } else if(length(indexes) == ncol(poset)) {
    # size1 + 1
    return(logsumexp(c(log_size1, log(1)))$logR )
  }
  
  log_size2 = genoLatticeSizeRecursion(as.matrix(poset[-indexes,-indexes]) ) 
  logsumexp(c(log_size1, log_size2))$logR
}  





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

geno_compatible_fraction <- function(poset, obs_events, weights){
  compatible_geno = compatible_genotypes(obs_events, poset)
  
  C = sum(weights[compatible_geno$compatible_indexes])
  N = sum(weights)
  # I: number of incompatible genotypes
  I = N - C
  
  # fraction of the data that are compatible with the given poset
  C / N
}

candidate_posets <- function (obs_events, obs_weights, min_compatible_geno_fraction) {
  posets <- list()
  violations = violation_freqs_w(obs_events, obs_weights)
  epsilons = sort(unique(c(0, violations)))
  for (eps in epsilons) {
    poset = maximal_poset(violations, eps)
    
    if(geno_compatible_fraction(poset, obs_events, obs_weights)  < min_compatible_geno_fraction) 
      break;
    
    colnames(poset) = rownames(poset) = colnames(obs_events)
    if (is_identical(poset, posets) == FALSE) {
      posets[[length(posets) + 1]] = poset
    }
  }
  posets
}




compute_q_e <- function(poset) {
  p=ncol(poset)
  if(ncol(poset) < 40 || sum(poset) < 8) {
    log_lattice_size = genoLatticeSize_fast(poset)
    log_q_e = -logsumexp( c(p * log(2), log_lattice_size), c(1, -1) )$logR
  } else {
    log_q_e = -  ncol(poset) * log(2)
  }
  #max(log_q_e, -  ncol(poset)/3 * log(2) )
  log_q_e
}


incompatible_loglike <- function(poset, obs_events, sampling_times, weights, compatible_geno, noise_model, geno_prob_noise=NULL) {
  p = ncol(poset)
  
  # C: number of compatible genotypes
  C = sum(weights[compatible_geno$compatible_indexes])
  N = sum(weights)
  # I: number of incompatible genotypes
  I = N - C
  
  # fraction of the data that are compatible with the given poset
  alpha = C / N
  
  if(noise_model == "uniform") {
    log_q_e = I * compute_q_e(poset)
    
  } else if(noise_model == "empty") {
    
    nc_indexes = setdiff(1:nrow(obs_events), compatible_geno$compatible_indexes)
    
    log_q_e = 0
    
    if(length(nc_indexes) > 0) {
      log_q_e = sum( weights[nc_indexes] * geno_prob_noise)
    } 
    
  } else {
    stop("Unknown noise model!")
  }
  
  incompatible_ll = 0.0
  if( I > 0 ){
    incompatible_ll = I * log(1-alpha)  +  log_q_e
  }
  
  list(ll=incompatible_ll,  alpha=alpha) 
}
