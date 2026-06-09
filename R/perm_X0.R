#' Permutation of X0
#'
#' Permutation of input matrix X0
#' @param X0 matrix; a score matrix X0, in which each column (layer) is a score vector over all vertices of A.
#' @param perms list of permutations of the rownames of X0, obtained by means of `permute_vertices()`
#' @return list with permutations of X0 (where the first element of the list contains real data)
#'
#' @export
#'

perm_X0 <- function(X0=NULL, perms=NULL){
	
	cat("perm_X0() is deprecated and will be replaced by perm_mat()\n")
	
	X0_list <- vector("list", length = length(perms))
	X0_list <- lapply(X0_list, function(x) X0)
	
	for(i in 2:length(X0_list)){ #assign the new rownames
		rownames(X0_list[[i]]) <- perms[[i]]
	}
	
	#reorder the rownames of the permutated X_i as the original matrix
	X0_list <- lapply(X0_list, function(X0_i) X0_i[match(rownames(X0_list[[1]]), rownames(X0_i)), , drop = F])
	
	return(X0_list)
	
}
