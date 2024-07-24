#' Permutation of graph vertices
#'
#' @param vert_deg named vector of vertex degree
#' @param k number of permutations.
#' @param seed_n optional numeric; to specify the seed for (pseudo) random number generation; using the same seed
#' @param method method to perform permutation; "simple" permutes vertices only based on vertex sets; "degree" permutes vertices based on vertex_sets and vertex degree
#' @param vertex_sets optional list with subsets of names(vert_deg); each subset will be permuteted separately; this is useful when distinct sets of vertices have different meaning;
#' @param cut_par number of intervals or width of the intervals
#' @param bin_type "interval" will use cut_interval(), "number" will use cut_number() while "width" will use cut_width() function from ggplot2
#' @return list with k+1 vectors of permutations of names(vert_deg), where the first element of the list contains original names(vert_deg)
#'
#' @export
#' @importFrom ggplot2 cut_interval cut_number cut_width
#' @importFrom stats setNames
#'

perm_vertices <- function(vert_deg=NULL, k=99, seed_n = NULL, vertex_sets=NULL, method=c("simple", "degree"), cut_par=20, bin_type=c("interval", "number", "width")){
  
  if(!is.null(seed_n)){
    set.seed(seed_n)
  }
  
  stopifnot(!is.null(names(vert_deg)))
  
  method <- match.arg(method, c("simple", "degree"))
  bin_type <- match.arg(bin_type, c("interval", "number", "width"))
  
  cat("Permutation type:", method , "\n")
  cat("Total permutations:", k+1, "\n")
  cat("Minimum possible FDR:", 1/(k+1), "\n")
  
  
  ### try to write the same code for all the possibilities, which therefore will depend on how objects are initialized, e.g. no bins means 1 bin
  ### degree preserving
  #assign every gene to a bin and pick up a representative from that bin
  
  bins <- rep(1, length(vert_deg))
  
  if(method == "degree"){
    
    cat("cut_par:", cut_par, "\n")
    cat("bin_type:", bin_type, "\n")

    if(bin_type == "interval"){
      bins <- as.numeric(cut_interval(vert_deg, cut_par)) ### use ggplot2
    }
    if(bin_type == "number"){
      bins <- as.numeric(cut_number(vert_deg, cut_par)) ### use ggplot2
    }
    if(bin_type == "width"){
      bins <- as.numeric(cut_width(vert_deg, cut_par)) ### use ggplot2
    }
  }
  
  if(is.null(vertex_sets)){
    
    vertex_sets <- setNames(rep(1, length(vert_deg)), names(vert_deg))
    
  }else{
    
    #remove shared vertices
    if(length(vertex_sets) > 1){
      
      for(i in 2:length(vertex_sets)){
        idx_shared <- which(vertex_sets[[i]] %in% vertex_sets[[i-1]])
        if(length(idx_shared)>0){
          vertex_sets[[i]] <- vertex_sets[[i]][-idx_shared]
        }
      }
      vertex_sets <-  vertex_sets[lengths(vertex_sets)>0]
      
      if(length(unique(unlist(vertex_sets))) != length(vert_deg)){
        stop("Vertex sets must cover all names of vert_deg.\n")
      }
    }
    
    vertex_sets <- setNames(rep(1:length(vertex_sets), times=lengths(vertex_sets)), unlist(vertex_sets))
    vertex_sets <- vertex_sets[match(names(vert_deg), names(vertex_sets))]
    stopifnot(length(vertex_sets) == length(vert_deg))
    
  }
  
  vertex_sets <- factor(paste(vertex_sets, bins, sep = "_"))
  cat("vertex sets and bins:", table(vertex_sets), "\n")
  vertex_sets_lev <- levels(vertex_sets)
  
  vertex_sets_list <- split(names(vert_deg), vertex_sets)
  
  l_fac <- unlist(lapply(lengths(vertex_sets_list), function(x) factorial(x)))
  cat("Min possible k:", min(l_fac), "\n")
  if(any(l_fac < k)){
    cat("Possible permutations: ", l_fac, "\n")
    cat("k:", k, "\n")
    stop("Possible permutations are less than requested, try reducing k, change cut_par or bin_type.\n")
  }
  
  ans <- vector("list", length = k+1)
  ans <- lapply(ans, function(x) names(vert_deg))
  
  for(i in 2:length(ans)){
    
    perm_i <- ans[[1]] #real permutation
    
    for(j in 1:length(vertex_sets_lev)){
      
      idx <- which(vertex_sets == vertex_sets_lev[j]) #rows of ith-level
      perm_i[idx] <- sample(perm_i[idx], length(perm_i[idx])) #permutaation of idx rows
      
    }
    
    ans[[i]] <- perm_i
  }
  
  return(ans)
  
}
