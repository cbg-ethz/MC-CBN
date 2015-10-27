

build_transition_matrix <- function( lambda, G) {
  .Call("build_transition_matrix", lambda, G)
}

probs <- function(lambda, G, times) 
{
  tr_matrix = t(build_transition_matrix(lambda, G) )
  
  
  indexes = which(tr_matrix != 0, arr.ind = T)
  tr_matrix <- new("dtTMatrix", x= tr_matrix[indexes], i= as.integer(indexes[,1]-1), 
                   j=as.integer(indexes[,2]-1), uplo='L', Dim= as.integer(c(nrow(tr_matrix),nrow(tr_matrix))))
  
  if(nrow(tr_matrix) > 100) {
    print(nrow(tr_matrix) )
  }
  v <- rep(0, nrow(tr_matrix))
  v[1] = 1

  t( sapply(times, function(x){ 
    A = expAtv(tr_matrix, v=v, t=x)$eAtv
    A[A<0]=0
    A/sum(A)
  }) )
}

genotype_probability <- function(poset, lambda, genotype, time) {
  G = compatible_genotypes_with_matrix(poset)
  
  geno_index = which(apply(G, 1, function(x) { all(genotype==x) }) )
  if(length(geno_index) == 0) {
    print(paste("The genotype (", paste(which(genotype==1), collapse=',') ,") is not compatible with the poset. Hence, the probability is zero") )
    return(0.0)
  }
  tr_matrix = build_transition_matrix(lambda, G) 
  expm(tr_matrix*time)[1, geno_index] 
}

mopt <- function(init_lambda, likelihood_fn, obs_events, dtimes, G, optim_method, p, control, init_eps, with_eps, ...) {
  res = list()
  lower = rep(1/mean(dtimes$all_times) * 10^(-5), p)
  upper = rep(Inf, p)
  init_params = init_lambda
  
  if(any(init_lambda < lower) ) {
    print("Initialized rates are not in the acceptable range! Another initialization method will be chosen!")
    init_lambda = rep(1/mean(dtimes$all_times), p)
  }
#   
#   if(with_eps==TRUE) {
#     init_params = c(init_lambda, init_eps)
#     lower = c(lower, 0.0)
#     upper = c(upper, 0.5)
#     
#   } 
  
  if(optim_method == "minqa") {
#     if(is.na(maxfun)) {
#       ctl = list(iprint=iprint)  
#     } else {
#       ctl = list(iprint=iprint, maxfun=maxfun)  
#     }
    
    fit = bobyqa(init_params, likelihood_fn, lower = lower, upper = upper, 
                 control = control, obs_events=obs_events,
                 dtimes = dtimes, G=G, ll_control = list(negate=TRUE, eps=init_eps), ...)
    res$value = -fit$fval
    res$par = fit$par
    res$fit = fit
    
  } else{
    fit= optim(par=init_params, fn=likelihood_fn, gr=NULL, ll_control = list(negate=TRUE, eps=init_eps), ..., obs_events=obs_events,
             dtimes = dtimes, G=G, method = c("L-BFGS-B"),
             lower = lower, upper = upper,
             control = control, hessian = FALSE)
    
    res$value = -fit$value
    res$par = fit$par
    res$fit = fit
  }
  res
}

#' fit_params
#' @export 
fit_params <- function(poset, obs_events, times, D=NA, ilambda = NA, use_grad=FALSE, optim_method="minqa", control=list(),
  ## TODO : check arguments, check if var of mutations are not zer 
  ## stop if the poset has (a)? loop
  ## nrow(poset)  != ncol(poset) 
  ## ncol(obs_event) != ncol(poset) 
  ## nrow(obs_events) != length(times)
  
  err_model=c("no-error", "hamming"), eps=NA, ...) {
  
  err_model = match.arg(err_model)
  
  if(use_grad==TRUE) {
    stop("TODO: no gradient is implemented for the likelihood! Try with grad=TRUE")
  }
  else  {
    grad_fn = NULL  
  }
  
  if(is.na(D)) {
    stop("TODO: implement a default dtimes without any approximation!")
  } else {
    dtimes = discretize_times(times, D)  
  }
  
  
  if(err_model == "no-error") {
    
    fit = estimate_lambda(obs_events, dtimes, poset, likelihood_fn = likelihood, init_lambda = ilambda, init_eps=eps, grad=grad_fn, 
                optim_method=optim_method, control=control, with_eps=FALSE, ...) 
    
  } else if(err_model == "hamming") {
    
    if(is.na(eps)) {
      stop("Please provide epsilon for hamming error modeling")
    }
    fit = estimate_lambda(obs_events, dtimes, poset, likelihood_fn = likelihood_with_eps, init_lambda = ilambda, init_eps=eps, grad=grad_fn, 
                    optim_method=optim_method, control=control, with_eps=TRUE, ...) 
  } else {
    stop(paste(err_model,": invalid error modeling!" ) )
  }
  
   
}



fit_params_with_eps <- function(poset, obs_events, times, D=NA, ilambda = NA, 
                                use_grad=FALSE, optim_method="minqa",  control=list(), tol=0.02, lower=0.0, upper=0.4, ...) {
  
  obj <- function(eps, poset, obs_events, times, D, ilambda, use_grad, optim_method, control,  ...) { 
#     print(eps)
    fit = fit_params(poset, obs_events, times, D=D, ilambda = ilambda, use_grad=use_grad, optim_method="optim_method", control=control,
                     err_model="hamming", eps=eps, ...)
    fit$value
  }
  
  epsilon <- optimize(obj, c(lower, upper), ..., poset=poset, obs_events=obs_events, times=times, D=D,  ilambda=ilambda, use_grad=use_grad, optim_method=optim_method, 
                      control=control,  tol = tol, maximum = TRUE)$maximum
  
  fit = fit_params(poset, obs_events, times, D=D, ilambda = ilambda, use_grad=use_grad, optim_method="optim_method", control=control,
                   err_model="hamming", eps=epsilon, ...)
  fit$epsilon = epsilon
  fit
}



#' estimate_lambda
#' @export 
estimate_lambda <- function(obs_events, dtimes, poset_, likelihood_fn = likelihood, init_lambda = NA, init_eps=0.3, grad=NULL, 
                         optim_method="minqa", control=list(), with_eps=FALSE, verbose=FALSE, ...) {
  p = nrow(poset_)
  if(any(is.na(init_lambda)) ) {
    init_lambda = initialize_lambda( obs_events, dtimes$all_times, poset_, verbose)    
  }
  
  graphobj <- graph.adjacency(poset_, mode = "undirected")
  membership = clusters(graphobj)$membership
  ll = 0
  
  mle_lambda = rep(0, ncol(poset_))
  fits = list()
  for (i in 1:max(membership)) {
    indexes = which(membership == i)
    poset = poset_[indexes, indexes]
    G = NA
    geno_indexes = NA
#     if (length(indexes) > 1) {
      G = compatible_genotypes_with_matrix(poset)
      geno_indexes = compute_geno_indexes(obs_events[, indexes], poset)
#     }
    
    res = mopt(init_lambda[indexes], likelihood_fn, obs_events[, indexes], dtimes, G,  
             optim_method, length(indexes), geno_indexes = geno_indexes,  control, init_eps, with_eps, ...) 
    ll = ll + res$value
    mle_lambda[indexes] = res$par
    
    fits[[length(fits)+1]] = res
  }
  
  list(value=ll, par=mle_lambda, fits=fits)
}

initialize_lambda <- function ( obs_events, times, poset, verbose=FALSE) {
  p = nrow(poset)
  
  if(p == 1) {
     return( 1/mean(times) )
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
  avg_time = mean(times)
  ###rep(1/mean(times),p)
  sapply(theta, function(x) { -log(max(1-x, 0.00001)) /avg_time })
}





genotype_probability_fast <- function(poset, lambda, genotype, time, nrOfSamples=100, topo_path=NULL) {
  
  if(is.null(topo_path) ) {
    topo_path = my.topological.sort(poset) -1
  }
  
#   G = compatible_genotypes_with_matrix(poset)
  
#   geno_index = which(apply(G, 1, function(x) { all(genotype==x) }) )
  if(is_compatible(genotype, poset) == FALSE) {
    print(paste("The genotype (", paste(which(genotype==1), collapse=',') ,") is not compatible with the poset. Hence, the probability is zero") )
    return(0.0)
  }

  T = .Call("drawHiddenVarsSamples", nrOfSamples, genotype, time, lambda,  poset, topo_path)
  values = cbn_density_(T$T, lambda) - T$densities
  logp = (logsumexp(values)$logR - log(nrOfSamples))
  exp(logp)
}