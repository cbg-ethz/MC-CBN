

estimate_mutation_rates <- function(poset, genotypes, sampling_times, weights = NULL, max_iter=100,  zeta = 0.2, ilambda = NULL,
                          nrOfSamples = 5, verbose = TRUE, maxLambdaValue=10^6) {
  MCEM(poset=poset, obs_events=genotypes, sampling_times=sampling_times, max_iter=max_iter, weights=weights, alpha = zeta, 
       nrOfSamples = nrOfSamples, verbose = verbose, maxLambdaValue=maxLambdaValue)   
}

# TODO: give error message if O is not matrix


MCEM <- function(poset, obs_events, sampling_times, max_iter=100,  weights=NULL, alpha = 0.2, ilambda=NULL, nrOfSamples = 50,
                 verbose = TRUE, small_number=10^-7, maxLambdaValue) {
  
  if( is_all_genotypes_compatible_with_poset(poset, obs_events, weights) == FALSE){
    stop("Error in the function MCEM: Some genotypes with non-zero weights are incompatible with the poset!")
  }

  if(is.null(weights)) {
    weights = rep(1, nrow(obs_events))
  } 

  if(is.null(ilambda)) {
    ilambda = initialize_lambda(obs_events = obs_events, sampling_times = sampling_times, poset = poset, verbose=verbose)
    ilambda[ilambda==0] = small_number
  }
  topo_path = my.topological.sort(poset)-1
  .Call("MCEM", ilambda, poset, obs_events, sampling_times, max_iter, alpha, topo_path, weights, nrOfSamples, verbose, maxLambdaValue)
}

allTopoSorts <- function(poset) {
  .Call("allTopoSorts", poset)
}

is_all_genotypes_compatible_with_poset <- function(poset, obs_events, weights) {
  
  if(is.null(weights)) {
    weights = rep(1, nrow(obs_events))
  }
  
  for (i in 1:nrow(obs_events)) {
    if (is_compatible(obs_events[i, ], poset) == FALSE && weights[i] > 0) {
      return(FALSE)
    }
  }
  return(TRUE)
}


# RcppExport SEXP MCEM( SEXP ilambda_, SEXP poset_, SEXP O_, SEXP times_, SEXP max_iter_, SEXP alpha_, NumericVector topo_path_) 

#MCEM( NumericVector lambda, NumericMatrix poset, NumericMatrix O, NumericVector times, SEXP max_iter_, SEXP alpha_) {
#MCEM(poset, O, times, max_iter=10, ilambda=NULL, alpha = 0.05, verbose = TRUE)


learn_network_boot <- function(obs_events, sampling_times, B = 50, weights=NULL, max_iter=200, alpha = 1.0, nrOfSamplesForEStep=50,
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
    
    
    result = learn_network(obs_events2, times2, weights2, max_iter, alpha, nrOfSamplesForEStep,
                  nrOfSamplesForLL, noise_model, verbose, min_compatible_geno_fraction)   
    ml_index = which.max(result$loglik)
    mle_poset = result$posets[[ml_index]]
    
    est_poset = est_poset + mle_poset
    est_poset_trans = est_poset_trans + trans_closure(mle_poset)
  }
  list(poset= est_poset/B, poset_trans = est_poset_trans/B)
}



learn_network <- function(obs_events, sampling_times, weights=NULL, max_iter=100, zeta = 0.2, L=5,
                          nrOfSamplesForLL = 100,   noise_model="empty", verbose=FALSE, min_compatible_geno_fraction=0.5, maxLambdaValue=10^6) {
  
  all_maximal_posets_mcem(obs_events, sampling_times, max_iter, zeta, L, nrOfSamplesForLL, weights, removeZeroWeights=TRUE,  
                          noise_model, verbose, min_compatible_geno_fraction, maxLambdaValue=maxLambdaValue)  
}


genotype_probs_empty_poset <- function (obs_events, sampling_times, weights, max_iter, alpha, nrOfSamplesForEStep, verbose  ) {
   # all the observed genotypes.
   poset_noise = make_empty_poset(ncol(obs_events) ) 
   
   # 1) estimate the rates from all the data use MCEM function
   fit = MCEM(poset_noise, obs_events, sampling_times, max_iter=max_iter,
              weights=weights, alpha = alpha, ilambda=NULL, nrOfSamples=nrOfSamplesForEStep, verbose=verbose)
   lambdas_noise = fit$par
   
   geno_prob_noise = all_genotype_prob_for_empty_poset(lambdas_noise, obs_events, sampling_times, log.p=TRUE)
   geno_prob_noise 
}

all_maximal_posets_mcem <- function(obs_events, sampling_times, max_iter=200, alpha = 1.0, nrOfSamplesForEStep=50,
                                    nrOfSamplesForLL = 100, weights=NULL, removeZeroWeights=TRUE, noise_model="empty", verbose=FALSE, 
                                    min_compatible_geno_fraction, maxLambdaValue)
{
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
  
  if(noise_model == "empty_approx") {
    geno_prob_noise = genotype_probs_empty_poset(obs_events, sampling_times, weights, max_iter, alpha, nrOfSamplesForEStep, verbose  )
  }
    
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
                                                     max_iter, alpha, nrOfSamplesForEStep, verbose)
        
      }
      # learn the empty poset params using incompatible data
    }

    logkile_incompatible = incompatible_loglike(poset, obs_events, sampling_times, weights, compatible_geno, noise_model, geno_prob_noise)
    
    if( (i > 1 && (max_loglike >  logkile_incompatible) ) || ( length(compatible_geno$compatible_indexes) <= 0 ) ) {
      cur_loglik = logkile_incompatible
      if(verbose==TRUE) {
        print("No parameter estimation for this poset. Upper bound likelihood of this poset is less than current maximum likelihood poset!")        
      }
    } else{
      fit = MCEM(poset, obs_events[compatible_geno$compatible_indexes, , drop=FALSE ], sampling_times[compatible_geno$compatible_indexes], max_iter=max_iter,
                 weights=weights[compatible_geno$compatible_indexes], alpha = alpha, ilambda=ilambda, nrOfSamples=nrOfSamplesForEStep, verbose=verbose, maxLambdaValue=maxLambdaValue)
      
      lambdas = fit$par
      if(verbose) {
        print(fit)
        print(proc.time() - t1)
      }
      t1 = proc.time()
      control = list(ll_method="importance", nrOfSamples=nrOfSamplesForLL )
      cur_loglik = loglike_mixture_model(poset, lambdas, obs_events, sampling_times, weights, control, compatible_geno, logkile_incompatible) 
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



################ computing likelihood using importance sampling

cbn_density_ <- function(T, rates) {
  p = ncol(T)
  apply(T, 1, function(x) { sum(dexp(x, rates, log=TRUE))}  )
}

loglike_importance_sampling <- function(poset, lambda, obs_events, sampling_times, nrOfSamples, weights, with_eps=FALSE, eps=NA)  {
  if(with_eps) {
    stop(paste("not implemented yet for with_eps=TRUE"))
  }
  topo_path = my.topological.sort(poset)-1
  probs = c()
  approx_loglike = 0.0
  
  for(i in 1:nrow(obs_events) ) {
    T = .Call("drawHiddenVarsSamples", nrOfSamples, obs_events[i,], sampling_times[i], lambda,  poset, topo_path)
    values = cbn_density_(T$T, lambda) - T$densities
    probs = c(probs, weights[i] * (logsumexp(values)$logR - log(nrOfSamples)))
    approx_loglike = approx_loglike + weights[i] * (logsumexp(values)$logR - log(nrOfSamples))
  }
  list(approx_loglike = approx_loglike , probs=probs)
}


genotype_prob_for_empty_poset <- function(lambda, genotype, time, log.p = FALSE) {
  p = length(lambda)
  
  mu_probs = sapply(1:p, function(i) {pexp(time, lambda[i], lower.tail = genotype[i] == 1, log.p=TRUE) } )
  
  log_sum = sum(mu_probs)
  if(log.p == TRUE)
    return(log_sum)
  exp(log_sum)
}

all_genotype_prob_for_empty_poset <- function(lambda, obs_events, sampling_times, log.p=FALSE) {
  
  probs = sapply(1:nrow(obs_events), function(i) { 
    genotype_prob_for_empty_poset(lambda, obs_events[i, ], sampling_times[i], log.p = log.p)
  } )
  probs
}
