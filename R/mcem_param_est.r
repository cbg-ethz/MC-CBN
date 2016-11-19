
initialize_lambda <- function ( obs_events, average_sampling_times, poset, verbose=FALSE) {
  p = nrow(poset)
  
  if(p == 1) {
    return( 1/average_sampling_times )
  }
  
  theta <- rep(0, p)
  for(i in 1:p) {
    parents = which(poset[, i] == 1)
    
    if(length(parents) == 0) {
      theta[i] = sum(obs_events[,i]) / nrow(obs_events)
    } else { 
      
      #       if(length(parents) == 1) {
      #         indexes = which( obs_events[, parents] ==  1)
      #       } else {
      #         indexes = which( apply(obs_events[, parents, drop=FALSE], 1, function(x) { all(x==1)}) )
      #       }
      
      indexes = which( apply(obs_events[, parents, drop=FALSE], 1, function(x) { all(x==1)}) )
      
      if(length(indexes) == 0) {
        if(verbose==TRUE) { 
          print(paste("No enough evidence for the mutation ", i, sep=''))
        }
        theta[i] = 0.000001
      } else {
        theta[i] = sum(obs_events[indexes,i]) / length(indexes)
      }
    } 
  }
  sapply(theta, function(x) { -log(max(1-x, 0.00001)) /average_sampling_times })
}



estimate_mutation_rates <- function(poset, genotypes, sampling_times=NULL, weights = NULL, max_iter=100,  zeta = 0.2, ilambda = NULL,
                          nrOfSamples = 5, verbose = TRUE, maxLambdaValue=10^6, lambda_s=1.0) {
  sampling_times_available = is.null(sampling_times) == FALSE
  MCEM(poset=poset, obs_events=genotypes, sampling_times=sampling_times, max_iter=max_iter, weights=weights, zeta = zeta, 
       nrOfSamples = nrOfSamples, verbose = verbose,  maxLambdaValue=maxLambdaValue, lambda_s = lambda_s, sampling_times_available=sampling_times_available)   
}

MCEM <- function(poset, obs_events, sampling_times=NULL, max_iter=100,  weights=NULL, zeta = 0.2, ilambda=NULL, nrOfSamples = 50,
                 verbose = TRUE, small_number=10^-7, maxLambdaValue, lambda_s, sampling_times_available) {
  
  if(sampling_times_available == FALSE) {
    average_sampling_times = lambda_s
    sampling_times = rep(0, nrow(obs_events))
  } else {
    average_sampling_times = mean(sampling_times)
  }

  if( is_all_genotypes_compatible_with_poset(poset, obs_events, weights) == FALSE){
    stop("Error in the function MCEM: Some genotypes with non-zero weights are incompatible with the poset!")
  }

  if(is.null(weights)) {
    weights = rep(1, nrow(obs_events))
  } 

  if(is.null(ilambda)) {
    ilambda = initialize_lambda(obs_events = obs_events, average_sampling_times = average_sampling_times, poset = poset, verbose=verbose)
    ilambda[ilambda==0] = small_number
  }
  topo_path = my.topological.sort(poset)-1
  .Call("MCEM", ilambda, poset, obs_events, sampling_times, max_iter, zeta, topo_path, weights, nrOfSamples, verbose, 
        maxLambdaValue, lambda_s, sampling_times_available)
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


