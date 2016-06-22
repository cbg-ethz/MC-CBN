
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


mixedRadixGeneration <- function (d, m)
{
  N <- prod(m)
  T <- matrix(0, N, d)
  a <- matrix(0, 1, d + 1)
  m <- cbind(2, m)
  for (i in 1:N) {
    T[i, ] = a[2:(d + 1)]
    j <- d
    while (a[j + 1] == m[j + 1] - 1) {
      a[j + 1] <- 0
      j <- j - 1
    }
    a[j + 1] = a[j + 1] + 1
  }
  return(T)
}


hyperCube <- function (n)
{
  m <- matrix(2, nrow = 1, ncol = n)
  H <- mixedRadixGeneration(n, m)
  return(H)
}


# The source code is copied partly from icbn package

#' orderIdeals
#' @export
orderIdeals <- function (poset)
{
  H <- hyperCube(ncol(poset))
  N_H <- nrow(H)
  p <- ncol(H)
  G <- matrix(1, nrow = N_H, ncol = 1)
  for (j1 in 1:p) {
    for (j2 in 1:p) {
      if (poset[j1, j2]) {
        G <- as.numeric(G & (H[, j1] >= H[, j2]))
      }
    }
  }
  list(G=G, H=H, lattice_size=sum(G==1))
}



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

# 
# #' compatible_genotypes
# #' @export
# compatible_genotypes <- function (X, poset)
# {
#   # H <- hyperCube(p)
#   G <- orderIdeals(poset)$G
#   
#   compatible_indexes = c()
#   N <- nrow(X)
#   p <- ncol(X)
#   E <- p - (1:p)
#   N_H <- 2^p
#   D <- rep(0, N_H)
#   for (i in 1:N) {
#     idx = 1 + sum(X[i, ] * 2^E)
#     D[idx] = D[idx] + 1
#     
#     if(G[idx] == 1)
#       compatible_indexes <- c(compatible_indexes, i)
#   }
#   
#   list( fraction=sum(D[as.logical(G)])/sum(D), compatible_indexes=compatible_indexes)
# }

# TODO: change the name to distinguish between this and the compatible_genotypes function
compatible_genotypes_with_matrix <- function(poset) {
  poset = as.matrix(poset)
  ordIdeals = orderIdeals(poset)
  possible_indexes = which(ordIdeals$G == 1)
  as.matrix(ordIdeals$H[possible_indexes, ], nrow=length(possible_indexes))
}


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
  
    # using Turan theorem for the maximum number of cover relations in a DAG
#   if(nr_edges > (nr_mut*nr_mut/4) ){
#     stop(paste("nr_edges should be less than ", nr_mut*nr_mut/4) )
#   }
  
  # Initialisation
  all_posets = rep(0, nr_pos)
  pos = 1
  
  while (pos <= nr_pos){
    
    poset <- matrix(0, nr_muts[pos], nr_muts[pos])
    poset[upper.tri(poset)][sample(choose(nr_muts[pos], 2), nr_edges[pos])] <- 1
    

#     poset = matrix(0, nr_mut, nr_mut)
# #     nr_edges = round(runif(1, 0, max_edges))
#     i = 1
#     iter = 0
#     while (i < nr_edges & iter < (nr_edges * 5) ){
#       
#       first_coordinate = round(runif(1,1,nr_mut))
#       second_coordinate = round(runif(1,1,nr_mut))
#       if(poset[first_coordinate, second_coordinate] == 0) {
#         
#         poset[first_coordinate, second_coordinate] = 1
#         g <- graph.adjacency(poset)
#         
#         if(is.dag(g)==FALSE){
#           poset[first_coordinate, second_coordinate] = 0
#         } else {
#           poset = trans_closure(poset)  
#           i = sum(poset)
#         }
#         #       else {i = i + 1}
#         
#         iter = iter + 1
#       }
#     }
      # added by Hesam. Add the transitive reduced poset to the list
    if(trans_reduced) {
      poset = trans_reduction(poset)  
    }
    
    all_posets[pos] = list(poset)
    
#     # Testing for uniqueness of new poset
#     test_uniqueness = rep(0, pos)
#     for (i in 1:pos){
#       test_uniqueness[i] = identical(all_posets[[pos]], all_posets[[i]])
#     }
#     If only identical with itself, poset is kept
#     if (sum(test_uniqueness)==1){
      pos = pos + 1  
#     }
  }
  all_posets
}

# 
# trans_reduction <- function(A) {
# #   R <- as.relation(A)
# #   
# #   RT = transitive_reduction(R)
# #   relation_incidence(RT)
#   hasse(A)
# }
# 
# 
# ############## Two following functions are copied from the ICBN package
# 
# hasse <- function (poset)
# {
#   p <- nrow(poset)
#   dimnames(poset) <- NULL
#   for (i in 1:p) {
#     for (j in 1:p) {
#       if (poset[i, j] == 1) {
#         poset[i, j] <- Inf
#         g <- graph(convert2EdgeList(poset), n = p)
#         dist <- shortest.paths(g, v = i, mode = "out",
#                                algo = "dijkstra")
#         if (is.finite(dist)[j] & (dist > 1)[j]) {
#           poset[i, j] <- 0
#         }
#         else {
#           poset[i, j] <- 1
#         }
#       }
#     }
#   }
#   poset <- apply(poset, 2, function(x) as.numeric(x > 0))
#   rownames(poset) <- colnames(poset)
#   return(poset)
# }



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

########### transitive closure and reduction of a poset



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



plot_poset <- function(robust_poset, size=12) {
  Names = colnames(robust_poset)
  colnames(robust_poset) <- Names
  rownames(robust_poset) <- Names
  am.graph <- new("graphAM", adjMat=robust_poset, edgemode="directed")
  plot(am.graph, attrs = list( node = list(color = "transparent", fontsize = size, fontcolor="dodgerblue4"), 
                               edge = list(arrowsize=0.5, color="antiquewhite4")))
}

# 
# logsumexp <- function(logs) {
#   m = max(logs)
#   S = sum(exp(logs-m))
#   m + log(abs(S))
# }


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
