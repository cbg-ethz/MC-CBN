#' @export
#' @noRd
my.topological.sort <- function(poset) {
  p = ncol(poset)
  sorted.list = rep(0, p)
  for(i in 1:p) {
    sorted.list[i] = min.degree = which.min(apply(poset, 2, sum))  
    poset[min.degree, ] = 0
    poset[, min.degree] = Inf
  }
  
  sorted.list
}

#' @export
#' @noRd
is_compatible <- function(genotype, poset) {
  p = nrow(poset)
  for(i in 1:p) {
    if(genotype[i] == 1) {
      parents = which(poset[, i] == 1)
      if(length(parents) > 0 ) {
        if(any(genotype[parents] == 0)) {
          return (FALSE)
        } 
      }
    }
  }
  return(TRUE)
}

# much faster than old implementation
#' @export
#' @noRd
compatible_genotypes <- function (X, poset)
{
  compatible_indexes = c()
  for(i in 1:nrow(X)) {
    if(is_compatible(X[i, ] , poset)) {
      compatible_indexes <- c(compatible_indexes, i)  
    }
  }
  list( fraction=length(compatible_indexes)/nrow(X), compatible_indexes=compatible_indexes)
}

#' @export
#' @noRd
random_poset <- function(p, graph_density=0.15, trans_reduced = TRUE){
  random_posets(nr_pos=1, nr_muts=p , ldenses=graph_density, trans_reduced = trans_reduced)[[1]]
}

random_posets <- function(nr_pos, nr_muts , ldenses, trans_reduced = TRUE){

  if(length(ldenses) == 1 ) {
    ldenses = rep(ldenses, nr_pos)
  } else if(length(ldenses) < nr_pos ) {
    stop( paste("Invalid value for ldenses ", ldenses) )
  }
  
  if(length(nr_muts) == 1 ) {
    nr_muts = rep(nr_muts, nr_pos)
  } else if(length(nr_muts) < nr_pos ) {
    stop( paste("Invalid value for nr_muts ", nr_edges) )
  }
  
  nr_edges = ceiling(choose(nr_muts, 2)*ldenses)
  
  # Initialisation
  all_posets = rep(0, nr_pos)
  pos = 1
  
  while (pos <= nr_pos){
    
    poset <- matrix(0, nr_muts[pos], nr_muts[pos])
    poset[upper.tri(poset)][sample(choose(nr_muts[pos], 2), nr_edges[pos])] <- 1

      # added by Hesam. Add the transitive reduced poset to the list
    if(trans_reduced) {
      poset = trans_reduction(poset)  
    }
    
    all_posets[pos] = list(poset)
    pos = pos + 1  
  }
  all_posets
}

########### transitive closure and reduction of a poset
#' @export
#' @noRd
trans_closure <- function(A) {
  old_names = dimnames(A)
  
  colnames(A) =  rownames(A) = paste("n", 1:nrow(A), sep='')
  R <- as.relation(A)
  
  RT = transitive_closure(R)
  res = relation_incidence(RT)
  res = res[rownames(A),colnames(A)]
  
  dimnames(res) = old_names
  res
}

#' @export
#' @noRd
trans_reduction <- function(A) {
  
  old_names = dimnames(A)
  
  colnames(A) =  rownames(A) = paste("n", 1:nrow(A), sep='')
  R <- as.relation(A)
  
  RT = transitive_reduction(R)
  res = relation_incidence(RT)
  res = res[rownames(A),colnames(A)]
  
  dimnames(res) = old_names
  res
}


############# plotting a cbn
#' @export
#' @noRd
plot_poset <- function(robust_poset, size=12) {
  Names = colnames(robust_poset)
  colnames(robust_poset) <- Names
  rownames(robust_poset) <- Names
  am.graph <- new("graphAM", adjMat=robust_poset, edgemode="directed")
  plot(am.graph, attrs = list( node = list(color = "transparent", fontsize = size, fontcolor="dodgerblue4"), 
                               edge = list(arrowsize=0.5, color="antiquewhite4")))
}

##########


logsumexp <- function(logs, signs=NULL) {
  if(is.null(signs)) {
    signs = rep(1, length(logs) )
  }
  indexes = which(signs != 0)
  signs = signs[indexes]
  logs = logs[indexes]

  m = max(logs)
  S = sum(signs*exp(logs-m))
  logR = m + log(abs(S))
  list(logR=logR, sign=sign(S))
}

myExp <- function(myNum) {
  myNum$sign * exp(myNum$logR)
}

convert2EdgeList <- function(poset){
  p <- nrow(poset)
  edgeList <- list()
  if (is.null(rownames(poset))){
    for (i in 1:p){
      for (j in 1:p){
        if (poset[i,j]==1){
          edgeList <- c(edgeList, i, j)
        }
      }
    }
  } else {
    for (i in 1:p){
      for (j in 1:p){
        if (poset[i,j]==1){
          edgeList <- c(edgeList, rownames(poset)[i], rownames(poset)[j])
        }
      }
    }
  }
  edgeList <- as.vector(unlist(edgeList))
  return(edgeList)
}