#' Permutation of X0
#'
#' Permutation of input matrix X0
#' @param X0 matrix; a score matrix X0, in which each column (layer) is a score vector over all vertices of A.
#' @param A adjacency matrix; it could be either normalized or not
#' @param k number of permutations.
#' @param seed_n optional numeric; to specify the seed for (pseudo) random number generation; using the same seed
#' @param method method to perform permutation of the row names of X0; "simple" permutes row names only based on vertex sets; "degree" permutes row names based on vertex_sets and vertex degree
#' @param vertex_sets NULL; list with subsets of X0 rownames; to control the permutation of values of each column of X0 within particular sets of vertices; this is useful when the values of at least one column of X0 can be assigned only to a subset of vertices;
#' @param cut_par number of intervals or width of the intervals
#' @param bin_type "interval" will use cut_interval(), "number" will use cut_number() while "width" will use cut_width() function from ggplot2
#' @return list with permutations of X0 (where the first element of the list contains real data)
#'
#' @export
#' @importFrom Matrix Matrix
#' @importFrom ggplot2 cut_interval cut_number cut_width
#' @importFrom stats setNames
#'

perm_X0 <- function(X0=NULL, A=NULL, k=99, seed_n = NULL, vertex_sets=NULL, method=c("simple", "degree"), cut_par=20, bin_type=c("interval", "number", "width")){
	
	if(!is.null(seed_n)){
		set.seed(seed_n)
	}
	
	method <- match.arg(method, c("simple", "degree"))
	bin_type <- match.arg(bin_type, c("interval", "number", "width"))
	
	stopifnot(rownames(X0)==rownames(A))
	
	cat("Permutation type:", method , "\n")
	cat("Total permutations:", k+1, "\n")
	cat("Minimum possible FDR:", 1/(k+1), "\n")
	
	
	X0 <- Matrix(X0)
	
	### try to write the same code for all the possibilities, which therefore will depend on how objects are initialized, e.g. no bins means 1 bin
	### degree preserving
	#assign every gene to a bin and pick up a representative from that bin
	
	row_bins <- rep(1, nrow(A))
	
	if(method == "degree"){
		cat("cut_par:", cut_par, "\n")
		cat("bin_type:", bin_type, "\n")
		d <- rowSums(sign(A))
		if(bin_type == "interval"){
			row_bins <- as.numeric(cut_interval(d, cut_par)) ### use ggplot2
		}
		if(bin_type == "number"){
			row_bins <- as.numeric(cut_number(d, cut_par)) ### use ggplot2
		}
		if(bin_type == "width"){
			row_bins <- as.numeric(cut_width(d, cut_par)) ### use ggplot2
		}
	}
	
	if(is.null(vertex_sets)){
		
		vertex_sets <- setNames(rep(1, nrow(A)), rownames(A))
		
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
			
			if(length(unique(unlist(vertex_sets))) != nrow(X0)){
				stop("Vertex sets must cover all rownames of X0.\n")
			}
		}
		
		vertex_sets <- setNames(rep(1:length(vertex_sets), times=lengths(vertex_sets)), unlist(vertex_sets))
		vertex_sets <- vertex_sets[match(rownames(A), names(vertex_sets))]
		stopifnot(length(vertex_sets) == nrow(A))
		
	}
	
	vertex_sets <- factor(paste(vertex_sets, row_bins, sep = "_"))
	cat("vertex sets and bins:", table(vertex_sets), "\n")
	vertex_sets_lev <- levels(vertex_sets)
	
	vertex_sets_list <- split(rownames(A), vertex_sets)
	
	l_fac <- unlist(lapply(lengths(vertex_sets_list), function(x) factorial(x)))
	if(any(l_fac < k)){
		cat("Possible permutations: ", l_fac, "\n")
		cat("k:", k, "\n")
		stop("Possible permutations are less than requested, try reducing k, change cut_par or bin_type.\n")
	}
	
	X0_list <- vector("list", length = k+1)
	X0_list <- lapply(X0_list, function(x) X0)
	
	for(i in 2:length(X0_list)){
		
		perm_i <- rownames(X0_list[[1]]) #real permutation
		
		for(j in 1:length(vertex_sets_lev)){
			
			idx <- which(vertex_sets == vertex_sets_lev[j]) #rows of ith-level
			perm_i[idx] <- sample(perm_i[idx], length(perm_i[idx])) #permutaation of idx rows
			
		}
		
		rownames(X0_list[[i]]) <- perm_i
	}
	
	
	X0_list <- lapply(X0_list, function(x) x[match(rownames(X0), rownames(x)), , drop = F])
	
	return(X0_list)
	
}
