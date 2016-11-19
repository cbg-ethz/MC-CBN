# 
# #' discretize_times
# #' @export 
# discretize_times <- function(times_, D) {
#   N = length(times_)
#   if(is.na(D)) {
#     D = ( max(times_) - min(times_) ) /50  
#   }
#   
#   t_indexes = rep(0, N)
#   
#   times2 = D * round(times_/D)
#   quan_times = unique(times2)
#   for(i in 1:N) {
#     t_indexes[which(times2 == quan_times[i])] = i  
#   }
#   
#   list(times = quan_times, t_indexes=t_indexes, all_times = quan_times[t_indexes])
# }
# 
# all_maximal_posets <- function(obs_events, sampling_times, D=NA, par=TRUE, num_par_posets=10, optim_method="minqa", control=list())
# {
#   if(ncol(obs_events) > 19) {
#     stop("Error, the exact method is only applicable to posets with less than 20 mutations!")
#   }
#   posets = candidate_posets( obs_events )
#   
#   nr_posets = length(posets)
#   
#   # initialize the variables
#   alphas= logliks = rep(0, nr_posets)
#   lambdas_mat = matrix(0, nr_posets, ncol(obs_events))
#   fits = list()
#   
#   if(par) {
#     mcoptions <- list(cores=num_par_posets)  
#     results = foreach(i = 1:nr_posets, .options.multicore = mcoptions) %dopar%
#     {
#       poset = posets[[i]]
#       compatible_geno = compatible_genotypes(obs_events, poset)
#       print(poset)
#       t1 = proc.time()
#       
#       dtimes = discretize_times(sampling_times[compatible_geno$compatible_indexes], D)
#       fit = estimate_lambda(obs_events[compatible_geno$compatible_indexes,], 
#                             dtimes, poset,
#                             init_lambda = NA, grad = NULL, 
#                             likelihood_fn=likelihood, maxfun=maxfun, optim_method=optim_method)
#       
#       lambdas = fit$par
#       print(fit)
#       print(proc.time() - t1)
#       t1 = proc.time()
#       #       cur_loglik = loglikelihood_cbn(poset, lambdas, sampling_times, obs_events, D)
#       ll_control = list(ll_method = "exact", D=D)
#       cur_loglik = loglike_mixture_model(poset, lambdas, obs_events, sampling_times, ll_control) 
#       print(cur_loglik)
#       print(proc.time() - t1)
#       print(paste("*** done ", sep = ""))
#       
#       list(fit=fit, alpha=cur_loglik$alpha, ll=cur_loglik$ll, lambdas=lambdas)
#     }
#     # update the variables for the poset i
#     for(i in 1:nr_posets) {
#       fits[[i]] = results[[i]]$fit
#       alphas[i] = results[[i]]$alpha
#       logliks[i] = results[[i]]$ll
#       lambdas_mat[i, ] = results[[i]]$lambdas
#     }
#   } else {
#     print("Serial computation: argument num_par_posets is ignored!")
#     for(i in 1:nr_posets) # for debugging simple for has to be used
#     {
#       poset = posets[[i]]
#       compatible_geno = compatible_genotypes(obs_events, poset)
#       print(poset)
#       t1 = proc.time()
#       dtimes = discretize_times(sampling_times[compatible_geno$compatible_indexes], D)
#       
#       fit = estimate_lambda(obs_events[compatible_geno$compatible_indexes,], 
#                             dtimes, poset, likelihood_fn=likelihood, 
#                             init_lambda = NA, grad = NULL, 
#                             optim_method=optim_method,  control=control)
#       
#       lambdas = fit$par
#       print(fit)
#       print(proc.time() - t1)
#       t1 = proc.time()
#       ll_control = list(ll_method = "exact", D=D)
#       cur_loglik = loglike_mixture_model(poset, lambdas, obs_events, sampling_times, ll_control) 
#       print(cur_loglik)
#       print(proc.time() - t1)
#       print(paste("*** done ", sep = ""))
#       
#       # update the variables for the poset i
#       fits[[i]] = fit
#       alphas[i] = cur_loglik$alpha
#       logliks[i] = cur_loglik$ll
#       lambdas_mat[i, ] = lambdas
#     }
#   }
#   
#   list(posets = posets, alphas = alphas, fits = fits,
#        logliks = logliks, lambdas_mat = lambdas_mat, obs_events = obs_events,
#        sampling_times = sampling_times)
# }
# 
# 
# 
# 
# build_transition_matrix <- function( lambda, G) {
#   .Call("build_transition_matrix", lambda, G)
# }
# 
# probs <- function(lambda, G, sampling_times) 
# {
#   tr_matrix = t(build_transition_matrix(lambda, G) )
#   
#   
#   indexes = which(tr_matrix != 0, arr.ind = T)
#   tr_matrix <- new("dtTMatrix", x= tr_matrix[indexes], i= as.integer(indexes[,1]-1), 
#                    j=as.integer(indexes[,2]-1), uplo='L', Dim= as.integer(c(nrow(tr_matrix),nrow(tr_matrix))))
#   
#   if(nrow(tr_matrix) > 100) {
#     print(nrow(tr_matrix) )
#   }
#   v <- rep(0, nrow(tr_matrix))
#   v[1] = 1
#   
#   t( sapply(sampling_times, function(x){ 
#     A = expAtv(tr_matrix, v=v, t=x)$eAtv
#     A[A<0]=0
#     A/sum(A)
#   }) )
# }
# 
# genotype_probability <- function(poset, lambda, genotype, time) {
#   G = compatible_genotypes_with_matrix(poset)
#   
#   geno_index = which(apply(G, 1, function(x) { all(genotype==x) }) )
#   if(length(geno_index) == 0) {
#     print(paste("The genotype (", paste(which(genotype==1), collapse=',') ,") is not compatible with the poset. Hence, the probability is zero") )
#     return(0.0)
#   }
#   tr_matrix = build_transition_matrix(lambda, G) 
#   expm(tr_matrix*time)[1, geno_index] 
# }
# 
# mopt <- function(init_lambda, likelihood_fn, obs_events, dtimes, G, optim_method, p, control, init_eps, with_eps, ...) {
#   res = list()
#   lower = rep(1/mean(dtimes$all_times) * 10^(-5), p)
#   upper = rep(Inf, p)
#   init_params = init_lambda
#   
#   if(any(init_lambda < lower) ) {
#     print("Initialized rates are not in the acceptable range! Another initialization method will be chosen!")
#     init_lambda = rep(1/mean(dtimes$all_times), p)
#   }
#   
#   if(optim_method == "minqa") {
#     fit = bobyqa(init_params, likelihood_fn, lower = lower, upper = upper, 
#                  control = control, obs_events=obs_events,
#                  dtimes = dtimes, G=G, ll_control = list(negate=TRUE, eps=init_eps), ...)
#     res$value = -fit$fval
#     res$par = fit$par
#     res$fit = fit
#     
#   } else{
#     fit= optim(par=init_params, fn=likelihood_fn, gr=NULL, ll_control = list(negate=TRUE, eps=init_eps), ..., obs_events=obs_events,
#                dtimes = dtimes, G=G, method = c("L-BFGS-B"),
#                lower = lower, upper = upper,
#                control = control, hessian = FALSE)
#     
#     res$value = -fit$value
#     res$par = fit$par
#     res$fit = fit
#   }
#   res
# }
# 
# #' fit_params
# #' @export 
# fit_params <- function(poset, obs_events, sampling_times, D=NA, ilambda = NA, use_grad=FALSE, optim_method="minqa", control=list(),
#                        ## TODO : check arguments, check if var of mutations are not zero 
#                        ## stop if the poset has (a)? loop
#                        ## nrow(poset)  != ncol(poset) 
#                        ## ncol(obs_event) != ncol(poset) 
#                        ## nrow(obs_events) != length(times)
#                        
#                        err_model=c("no-error", "hamming"), eps=NA, ...) {
#   
#   err_model = match.arg(err_model)
#   
#   if(use_grad==TRUE) {
#     stop("TODO: no gradient is implemented for the likelihood! Try with use_grad=FALSE")
#   }
#   else  {
#     grad_fn = NULL  
#   }
#   
#   if(is.na(D)) {
#     stop("TODO: implement a default dtimes without any approximation!")
#   } else {
#     dtimes = discretize_times(sampling_times, D)  
#   }
#   
#   
#   if(err_model == "no-error") {
#     
#     fit = estimate_lambda(obs_events, dtimes, poset, likelihood_fn = likelihood, init_lambda = ilambda, init_eps=eps, grad=grad_fn, 
#                           optim_method=optim_method, control=control, with_eps=FALSE, ...) 
#     
#   } else if(err_model == "hamming") {
#     
#     if(is.na(eps)) {
#       stop("Please provide epsilon for hamming error modeling")
#     }
#     fit = estimate_lambda(obs_events, dtimes, poset, likelihood_fn = likelihood_with_eps, init_lambda = ilambda, init_eps=eps, grad=grad_fn, 
#                           optim_method=optim_method, control=control, with_eps=TRUE, ...) 
#   } else {
#     stop(paste(err_model,": invalid error modeling!" ) )
#   }
#   
#   
# }
# 
# 
# fit_params_with_eps <- function(poset, obs_events, sampling_times, D=NA, ilambda = NA, 
#                                 use_grad=FALSE, optim_method="minqa",  control=list(), tol=0.02, lower=0.0, upper=0.4, ...) {
#   
#   obj <- function(eps, poset, obs_events, sampling_times, D, ilambda, use_grad, optim_method, control,  ...) { 
#     #     print(eps)
#     fit = fit_params(poset, obs_events, sampling_times, D=D, ilambda = ilambda, use_grad=use_grad, optim_method="optim_method", control=control,
#                      err_model="hamming", eps=eps, ...)
#     fit$value
#   }
#   
#   epsilon <- optimize(obj, c(lower, upper), ..., poset=poset, obs_events=obs_events, sampling_times=sampling_times, D=D,  ilambda=ilambda, use_grad=use_grad, optim_method=optim_method, 
#                       control=control,  tol = tol, maximum = TRUE)$maximum
#   
#   fit = fit_params(poset, obs_events, sampling_times, D=D, ilambda = ilambda, use_grad=use_grad, optim_method="optim_method", control=control,
#                    err_model="hamming", eps=epsilon, ...)
#   fit$epsilon = epsilon
#   fit
# }
# 
# 
# 
# #' estimate_lambda
# #' @export 
# estimate_lambda <- function(obs_events, dtimes, poset_, likelihood_fn = likelihood, init_lambda = NA, init_eps=0.3, grad=NULL, 
#                             optim_method="minqa", control=list(), with_eps=FALSE, verbose=FALSE, ...) {
#   p = nrow(poset_)
#   if(any(is.na(init_lambda)) ) {
#     init_lambda = initialize_lambda( obs_events, dtimes$all_times, poset_, verbose)    
#   }
#   
#   graphobj <- graph.adjacency(poset_, mode = "undirected")
#   membership = clusters(graphobj)$membership
#   ll = 0
#   
#   mle_lambda = rep(0, ncol(poset_))
#   fits = list()
#   for (i in 1:max(membership)) {
#     indexes = which(membership == i)
#     poset = poset_[indexes, indexes]
#     G = NA
#     geno_indexes = NA
#     #     if (length(indexes) > 1) {
#     G = compatible_genotypes_with_matrix(poset)
#     geno_indexes = compute_geno_indexes(obs_events[, indexes], poset)
#     #     }
#     
#     res = mopt(init_lambda[indexes], likelihood_fn, obs_events[, indexes], dtimes, G,  
#                optim_method, length(indexes), geno_indexes = geno_indexes,  control, init_eps, with_eps, ...) 
#     ll = ll + res$value
#     mle_lambda[indexes] = res$par
#     
#     fits[[length(fits)+1]] = res
#   }
#   
#   list(value=ll, par=mle_lambda, fits=fits)
# }
# 
# ## Likelihood fucntions
# .log_observed_geno_probs_decomp <- function(poset, lambda, obs_events, dtimes, with_eps=FALSE, eps=NA) {
#   graphobj <- graph.adjacency(poset, mode="undirected")
#   membership = clusters(graphobj)$membership
#   
#   lprobs = rep(0, nrow(obs_events) )
#   for(i in 1:max(membership)) {
#     indexes = which(membership == i)
#     G = geno_indexes = NA
#     if(length(indexes) > 1) {
#       G = compatible_genotypes_with_matrix(poset[indexes, indexes])
#       geno_indexes = compute_geno_indexes(obs_events[, indexes], poset[indexes, indexes])  
#     } 
#     if(with_eps==TRUE) {
#       G = compatible_genotypes_with_matrix(poset[indexes, indexes])
#       lprobs = lprobs + log_observed_geno_probs_with_eps(lambda[indexes], eps, obs_events[, indexes], 
#                                                          dtimes, G)
#     } else{ 
#       lprobs = lprobs + log_observed_geno_probs(lambda[indexes], obs_events[, indexes], 
#                                                 dtimes, G, geno_indexes)
#     }
#     
#   }
#   lprobs
# }
# 
# .likelihood_decomp <- function(poset, lambda, obs_events, dtimes) {
#   sum(.log_observed_geno_probs_decomp(poset, lambda, obs_events, dtimes))  
# }
# 
# loglike <- function(poset, lambda, obs_events, times, D, with_eps=FALSE, eps=NA)  {
#   dtimes = discretize_times(times, D)
#   
#   sum(.log_observed_geno_probs_decomp(poset, lambda, obs_events, dtimes, with_eps, eps) )
# }
# 
# 
# compute_geno_indexes <- function(obs_events, poset) {
#   if(is.matrix(obs_events) == FALSE) {
#     obs_events = as.matrix(obs_events)
#   }
#   
#   
#   G = compatible_genotypes_with_matrix(poset)
#   
#   geno_indexes = rep(0, nrow(obs_events))
#   for(i in 1:nrow(obs_events)) {
#     genotype = obs_events[i, ]
#     geno_index = which(apply(G, 1, function(x) { all(genotype==x)} ))
#     geno_indexes[i] = ifelse(length(geno_index) == 0, NA, geno_index)
#   }
#   geno_indexes
# }
# 
# 
# log_observed_geno_probs <- function(lambdas, obs_events, dtimes, G, geno_indexes) {
#   sampling_times = dtimes$times
#   t_indexes = dtimes$t_indexes
#   all_times = dtimes$all_times
#   if(length(lambdas) == 1) {
#     log_observed_geno_probs = log( abs(obs_events - exp(-lambdas*all_times) ) )
#     log_observed_geno_probs[is.infinite(log_observed_geno_probs)] = -50
#     return(  log_observed_geno_probs )
#   }
#   
#   if(any(is.na(geno_indexes) ) ) {
#     
#     print("Some genotypes are not compatible with the poset. Hence, the log likelihood is a very small number (minus infinity)")
#     # -Inf does not work with numerical optimization. Hence, -50 multiply by number of genotypes is used for very small likelihood 
#     return(rep(-50, length(all_times)))
#     
#   }
#   
#   # computing the genotype probabilities for all time points (and all possible genotypes). The return value
#   # is a matrix 
#   geno_probs = probs(lambdas, G, sampling_times)
#   
#   # computing the genotype probabilities for all time points only for  the observed genotypes. The return value
#   # is a vector
#   log_observed_geno_probs = log( sapply(1:nrow(obs_events), function(i) { geno_probs[t_indexes[i], geno_indexes[i]] }) )
#   
#   #   print(log_observed_geno_probs)
#   # Handling a special case: if the probability for a compatible genotype is zero (for any reason)
#   log_observed_geno_probs[is.infinite(log_observed_geno_probs)] = -50
#   log_observed_geno_probs   
# }
# 
# 
# likelihood <- function(lambdas, obs_events, dtimes, G, geno_indexes, ll_control=list(negate=FALSE)) {
#   negate = ll_control$negate
#   C = ifelse(negate, -1, 1)
#   
#   C*sum(log_observed_geno_probs(lambdas, obs_events, dtimes, G, geno_indexes))  
# }
# 
# 
# #####################################  add error modeling
# log_observed_geno_probs_with_eps<- function(lambda, eps, obs_events, dtimes, G) {
#   obs_events = as.matrix(obs_events, ncol=length(lambda))
#   sampling_times = dtimes$times
#   t_indexes = dtimes$t_indexes
#   all_times = dtimes$all_times
#   
#   # computing the genotype probabilities for all time points (and all possible genotypes). The return value
#   # is a matrix 
#   geno_probs = probs(lambda, G, sampling_times)
#   #   range(geno_probs)
#   prob_Y = rep(0, nrow(obs_events))
#   for( i in 1:nrow(obs_events)) {
#     Y = obs_events[i, ]
#     
#     P_Y_X = apply(G, 1, Y=Y, eps=eps, function(x, Y, eps) { 
#       prob_hamming_distance(x, Y, eps)
#     })
#     
#     prob_Y[i] = sum(P_Y_X * geno_probs[t_indexes[i], ] )
#   }
#   
#   
#   log_observed_geno_probs = log(prob_Y)
#   
#   # Handling a special case: if the probability for a genotype is zero (for any reason)
#   log_observed_geno_probs[is.infinite(log_observed_geno_probs)] = -50
#   
#   log_observed_geno_probs
# }
# 
# likelihood_with_eps <- function(lambda, obs_events, dtimes, G, geno_indexes, ll_control=list(negate=FALSE, eps=0.25)) {
#   
#   negate = ll_control$negate
#   eps = ll_control$eps
#   
#   log_observed_geno_probs = log_observed_geno_probs_with_eps(lambda, eps, obs_events, dtimes, G)
#   C = ifelse(negate, -1, 1)
#   C * sum(log_observed_geno_probs)
# }
# 
# 
# 
# 
# 
# # Y: observed genotype
# # X: true genotype
# # eps: per-locus error rate
# prob_hamming_distance <- function(X, Y, eps) {
#   p = length(X)
#   d = sum( X != Y )
#   (eps^d) * (1-eps)^(p-d)
# }



#########################  sampling time prediction
# .tlikelihood <- function(t, lambdas, m) {
#   P = pexp(t, lambdas)
#   
#   logP = m*P + (1-m)*(1-P)
#   logP[ logP==0 ] = 10^-50
#   -sum( log( logP ) )
# }
# 
# .fit_exp <- function(x, lambda_s=1, verbose=TRUE) {
#   M = sum(x == 1)
#   N = length(x)
#   
#   if(M == 0 || N==M) {
#     stop("Error! M == 0 or N==M. Get out of my sight!")
#   }
#   
#   lambda = M * lambda_s / (N-M)
#   expected = N* (lambda/(lambda+lambda_s) )
#   list(lambda=lambda,expected = expected)
# }
# 
# 
# .MLE_sampling_time_for_genotype <- function(x, rates, itime=1) {
#   optim(itime, .tlikelihood,  method = "L-BFGS-B", lower = 0.000001, upper=50, lambdas=rates, m=x)$par
# }
# 
# estimate_sampling_times <- function(X, N) {
#   # we assume a naive Bayes model. We first esimate the rate of each mutation.
#   rates = apply(X, 2, function(x) { .fit_exp(x, lambda_s=1, verbose=TRUE)$lambda}) 
#   
#   # Then we find the ML estimation of sampling time for each genotype with the exponential rates
#   ranges = t(apply(X, 1, expected_sampling_time_interval_for_genotype, rates=rates, N))
#   
#   apply(ranges, 1, function(x) {
#     if(is.infinite(x[2]) ) {
#       return(1.1*x[1] )
#     } else{
#       return(mean(x))
#     }
#   })
# }
# 
# expected_sampling_time_interval_for_genotype <- function(geno_, rates, N=1000) {
#   geno = which(geno_ == 1)
#   T = c()
#   for(rate in rates) {
#     T = cbind(T, rexp(N, rate) )
#   }
#   
#   if(length(geno) == 0) {
#     T_obs = 0
#     T_unobs = mean(apply(T, 1, min)  )  
#   } else {
#     T_obs = mean(apply(T, 1, function(x){ max(x[geno]) }) )  
#     
#     if(length(geno) != length(rates)) {
#       T_unobs = mean(apply(T, 1, function(x){ min(x[-geno]) }) )  
#     } else  {
#       T_unobs = Inf
#     }
#     
#   }
#   #     
#   #   print(T_obs)
#   #   print(T_unobs)
#   #   print(T_obs + T_unobs)
#   #   
#   #   print( )
#   c(T_obs, T_obs + T_unobs)
# }
# 
# 
# estimate_sampling_times2 <- function(X) {
#   # we assume a naive Bayes model. We first esimate the rate of each mutation.
#   rates = apply(X, 2, function(x) { .fit_exp(x, lambda_s=1, verbose=TRUE)$lambda}) 
#   
#   # Then we find the ML estimation of sampling time for each genotype with the exponential rates
#   apply(X, 1, .MLE_sampling_time_for_genotype, rates=rates, itime=1)
# }
# 


################ likelihood
# 
# #' loglike_mixture_model
# #' @export
# loglike_mixture_model <- function(poset, lambda, obs_events, sampling_times, weights, control = list(ll_method="importance"), compatible_geno, incomp_loglike){
#   
#   incompatible_ll = incomp_loglike$ll
#   C = sum(weights[compatible_geno$compatible_indexes])
#   alpha = incomp_loglike$alpha
#   
#   compatible_ll = 0.0
#   if( C > 0 )
#   {
#     genotypes = obs_events[compatible_geno$compatible_indexes, , drop=FALSE ]
#     
#     sampling_times = sampling_times[compatible_geno$compatible_indexes]
#     weights = weights[compatible_geno$compatible_indexes]
#     
#     if(control$ll_method == "importance") {
#       tmp = loglike_importance_sampling(poset, lambda, genotypes, sampling_times, control$nrOfSamples, weights, with_eps=FALSE, eps=NA)
#       compatible_ll = tmp$approx_loglike
#     } else {
#       dtimes = discretize_times(sampling_times, control$D)
#       compatible_ll = .likelihood_decomp(poset, lambda, genotypes, dtimes)
#     }
#     
#     compatible_ll = compatible_ll + C * log(alpha)
#   }
#   
#   list(ll=incompatible_ll + compatible_ll, incompatible_ll=incompatible_ll,  compatible_ll=compatible_ll, alpha=alpha)  
# }

#####

# 
# # TODO: change the name to distinguish between this and the compatible_genotypes function
# compatible_genotypes_with_matrix <- function(poset) {
#   poset = as.matrix(poset)
#   ordIdeals = orderIdeals(poset)
#   possible_indexes = which(ordIdeals$G == 1)
#   as.matrix(ordIdeals$H[possible_indexes, ], nrow=length(possible_indexes))
# }

################### from common

# 
# ########### others
# allTopoSorts <- function(poset) {
#   .Call("allTopoSorts", poset)
# }
# 
# 
# 
# 
# mixedRadixGeneration <- function (d, m)
# {
#   N <- prod(m)
#   T <- matrix(0, N, d)
#   a <- matrix(0, 1, d + 1)
#   m <- cbind(2, m)
#   for (i in 1:N) {
#     T[i, ] = a[2:(d + 1)]
#     j <- d
#     while (a[j + 1] == m[j + 1] - 1) {
#       a[j + 1] <- 0
#       j <- j - 1
#     }
#     a[j + 1] = a[j + 1] + 1
#   }
#   return(T)
# }
# 
# 
# hyperCube <- function (n)
# {
#   m <- matrix(2, nrow = 1, ncol = n)
#   H <- mixedRadixGeneration(n, m)
#   return(H)
# }
# 
# 
# # The source code is copied partly from icbn package
# 
# #' orderIdeals
# #' @export
# orderIdeals <- function (poset)
# {
#   H <- hyperCube(ncol(poset))
#   N_H <- nrow(H)
#   p <- ncol(H)
#   G <- matrix(1, nrow = N_H, ncol = 1)
#   for (j1 in 1:p) {
#     for (j2 in 1:p) {
#       if (poset[j1, j2]) {
#         G <- as.numeric(G & (H[, j1] >= H[, j2]))
#       }
#     }
#   }
#   list(G=G, H=H, lattice_size=sum(G==1))
# }
# 
# 
# 
# 

