#############################################################################################################
#############################################################################################################
############################################ independence Tests #############################################
#' weibull_loglike
#' @export
weibull_loglike <- function(x, E, T ) {
  x1 <- x[1]
  x2 <- x[2]
  
  ll = sum(log(pweibull(T[which(E == 1)], shape=x1, scale=x2) ) ) + sum(log(1-pweibull(T[which(E == 0)], shape=x1, scale=x2) ) ) 
  if(is.nan(ll)  || is.infinite(ll)) {
    return(length(E) * 100)
  }
  
  -ll
}

#' fit_weibull
#' B * N/4 is the maximum allwed threshold in the last iteration. The error threshold in iteration "iter" is defined
#' as max( 0.2, B * (N/4) * iter/max_iteration). The B * N/4 still gives a OK result. B * (N/4) is an approximation for
#' B * fa * fb/ (fa+fb) where fa and fb are the frequncy of events a and b. The derivation of this formula is straightforward.
#' 0.1 is a reasonable value for this parameter.
#' @export
fit_weibull <- function(E, T, org_init_par=NULL, B=0.1, verbose=F) {
  if(is.null(org_init_par)) {
    org_init_par = tryCatch( fitdistr(T, "weibull")$estimate,  error=function(e) c(mean(T)/1000, 10^6))  
#     org_init_par = c(mean(T), 1)
    
  }
  init_par = MLE_par =  org_init_par
  
  MLE_P = NA
  min_err = length(E) # maximum error for initialization
  
  i = 0
  N = length(E)
  max_iteration= 200
  
  # An error below the following threshold does not have a big impact on final pvalue
  error_thr = max(B * (N/4) * (i+1)/max_iteration, 0.2)
  
  # repeat until the global maximum likelihood solution is obtained or maximum number of iteration is reached.
  # The error is defined as the difference between expected and observed number of mutations. The error threshold 
  # slightly increases in each iteration.
  while( i < max_iteration & (min_err > error_thr ) ) {
    if(is.infinite( weibull_loglike(init_par, E, T) ) == FALSE) {
      #       opt1 = optim(init_par, fr, E=E, T=T, method="L-BFGS-B", lower=10^-9 * mean(stimes) )
      opt1 = optim(init_par, weibull_loglike, E=E, T=T)
      
      P = pweibull(T, shape=opt1$par[1], scale=opt1$par[2]) 
      err = abs(sum(P) -sum(E))
      
      if(err < min_err) {
        MLE_P = P
        min_err = err
        MLE_par = opt1$par
      }
      i = i + 1
      error_thr = max(B * (N/4) * (i+1)/max_iteration, 0.2)
    } 
    
    init_par[1] = runif(1, org_init_par[1]/50, org_init_par[1]*5)
    init_par[2] = runif(1, org_init_par[2]/10, org_init_par[2]*5)
  }
  
  if(verbose) {
    print(i)
    if(min_err > error_thr) {
      print("No convergence!")
      print(min_err)
      
    } 
    print("parameters")
    print(MLE_par)
    print(min_err)
  }
  
  MLE_P
}

expected_frequency <- function(X, Y, T, b, verbose) {
  Px = fit_weibull(X, T, verbose=verbose, B=b) 
  Py = fit_weibull(Y, T, verbose=verbose, B=b) 
  
  c(sum((1-Px) * (1-Py) ), sum(Px * (1-Py)),
    sum((1-Px) * Py), sum(Px * Py) )  
}


# X, Y  : mutation vectors
# T : sampling times
# S: conditional set
# B: see the description for fit_weibull
# adaptDF: to change the degree of freedom for empty conditional sets (zero counts)
CondGTest <- function(X, Y, T, S=NULL, B, adaptDF=FALSE, verbose=FALSE) {  
  if(is.null(S)) {
    S_indexes = rep(1, length(X))
  } else{
    if(ncol(S) > 1) {
      S_indexes = apply(S, 1, function(s) { sum(s* 2^{(length(s):1)-1}  ) } )  
    } else {
      S_indexes = S
    }
    
  }
  nzero_indexes = unique(S_indexes)
  stat_total = 0.0
  for(k in nzero_indexes ) {
    Xc = X[S_indexes == k]
    Yc = Y[S_indexes == k]
    Tc = T[S_indexes == k]
    Oc = c(sum(Xc==0 & Yc==0), sum(Xc==1 & Yc==0), sum(Xc==0 & Yc==1), sum(Xc==1 & Yc==1))
    
    Ec = expected_frequency(Xc, Yc, Tc, B, verbose)
    
    if(verbose) {
      print("expected")
      print(Ec)
      print("observed")
      print(Oc)
    }
    tmp = 2*Oc*log(Oc/Ec) 
    tmp[which(is.nan(tmp))] <- 0
    stat_total = stat_total + sum(tmp)  
  }
  if(is.null(S) ) {
    df=1
  } else {
    if(is.matrix(S)) {
      lenS <- ncol(S)    
    } else {
      lenS = 1
    }
    
    df <- 2^lenS
  }
  if(adaptDF)
    df <- length(nzero_indexes)
  if(verbose) {
    print(stat_total)  
  }
  
  pchisq(stat_total, df, lower.tail = FALSE)
}


CondGTest_pcform <- function(x, y, S, suffStat, B=0.1, adaptDF=TRUE, verbose=FALSE)  {
  
  if(length(S) == 0) {
    Sm= NULL
  } else {
    Sm = as.matrix(suffStat[, S])
  }
  
  CondGTest(suffStat[, x], suffStat[, y], suffStat[, "data.Ts"], Sm,  B=B, adaptDF = adaptDF, verbose=verbose )
}

#############################################################################################################
#############################################################################################################
#################################### Poset learning using conditional tests #################################


############################# find partial dag
# find partial dag using pc algorithm
find_pdag <- function(mutations, times, alpha, mmax=Inf, test =CondGTest_pcform, verbose=TRUE ) {
  
  n <- nrow(data)
  V <- colnames(mutations)
  if(is.null(V)) {
    V <- paste("V", 1:ncol(mutations), sep='')  
  }
  
  
  if(identical(test, disCItest) ) {
    suffStat <- list(dm = mutations, nlev =rep(2, ncol(mutations)), adaptDF = FALSE)
  } else  if(identical(test, binCItest) ) { 
    suffStat <- list(dm = mutations, adaptDF = FALSE)
  }  else{
    suffStat = data.frame(data=cbind(mutations, Ts=times))
  }
  
  fit <- pc(suffStat = suffStat,
            indepTest = test, ## indep.test: partial correlations
            alpha=alpha,  labels = V, m.max=mmax, verbose = verbose)
  
  fit
}



############################# orient the pdag
# a simple method for finding the temporal order among two mutations
# return value:
#   1 : X->Y
#   2 : Y->X
#   3 : independence
# 
orient_pair <- function(X, Y, T) {
  a = table(X, Y)[1,2]
  b = table(X, Y)[2,1]
  c = table(X, Y)[2,2]
  
  if(c == min(a, b, c)) {
    return(3) # independence
  }
  # less violations of the model X -> Y  
  if(a < b) {
    return(1)
  }
  return(2)
}


orient_pdag <- function(pdag, mutations, times) {
  p = ncol(mutations)
  fdag = matrix(0, p, p)
  for(i in 2:p) {
    for(j in 1:(i-1)) {
      if(max(pdag[i, j], pdag[i, j]) > 0) {
        or_ij = orient_pair(mutations[, i], mutations[, j], times)
        
        bcycle = FALSE
        if(or_ij == 1) {
          fdag[i, j] = 1
          if(is.dag(graph.adjacency(fdag))  == FALSE) {
            fdag[i, j] = 0
            bcycle = TRUE
          }
        } 
        
        if(or_ij == 2 || bcycle) {
          fdag[j, i] = 1
          if(is.dag(graph.adjacency(fdag))  == FALSE) {
            fdag[j, i] = 0
          }
        }
      }
      
    }
  }
  fdag
}

learn_poset_pc <- function(mutations, times, alpha, mmax=Inf, test =CondGTest_pcform, verbose=TRUE ) {
  if(is.null( names(mutations) ) ) { 
    colnames(mutations) = paste("M", 1:ncol(mutations), sep='')
  }
  pdag = find_pdag(mutations, times, alpha, mmax=mmax, test =CondGTest_pcform, verbose=TRUE ) 
  poset = convert_to_adja_matrix(pdag@graph)
  list(poset=trans_reduction(orient_pdag(poset, mutations, times)), pdag=pdag)
}



convert_to_adja_matrix <- function(graph) {
  gr = matrix(0, length( graph@nodes), length( graph@nodes))
  edges = matrix(0, 0 , 2)
  for(i in 1:length(graph@edgeL) ) {
    for(j in graph@edgeL[[i]]$edges) {
      edges = rbind(edges, c(i, j) )
      gr[i,j] = gr[j, i] = 1
    }
  }
  gr
}


learn_struct <- function( obs_events, times, method=c("pc", "pc-bagging"), ....) {
  if(method=="pc") {
    res = learn_poset_pc(obs_events, times, ... )
  } else if(method=="pc-bagging") {
    res = learn_poset_pc_bagging(obs_events, times, ... )
  }
  
  res
}



learn_poset_pc_bagging <- function(obs_events_, times_, B, thr=0.5, alpha=0.05, mmax=1, test =CondGTest_pcform, verbose=TRUE ) {
  aggregated_poset = matrix(0, ncol(obs_events_), ncol(obs_events_))
  
  for(i in 1:B) {
    print(paste("################ i=", i ,"#########################") )
    
    indexes = sample(nrow(obs_events_),replace=TRUE)
    obs_events = obs_events_[indexes, ]
    times = times_[indexes]
    
    if(is.null( names(obs_events) ) ) { 
      colnames(obs_events) = paste("M", 1:ncol(obs_events), sep='')
    }
    pdag = find_pdag(obs_events, times, alpha, mmax=mmax, test =CondGTest_pcform, verbose=TRUE ) 
    poset = convert_to_adja_matrix(pdag@graph)
    poset = orient_pdag(poset, obs_events, times) 
    
#      poset = trans_reduction(poset)
#      poset = trans_closure(poset)
#     poset = poset + t(poset)
     poset[poset>1]=1

    aggregated_poset = aggregated_poset + poset
  }
#   poset[poset < thr] = 0
#   poset[poset >= thr] = 1
  
  #trans_reduction() 
  aggregated_poset/B
}
 