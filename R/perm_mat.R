#' Permutation of X
#'
#' Permutation of input matrix X
#' @param X matrix; a score matrix X, in which each column (layer) is a score vector over all vertices of A.
#' @param perms list of permutations of the rownames of X, obtained by means of `permute_vertices()`.
#' @param columns whether to apply the permutation to the columns.
#' @return list with permutations of X (where the first element of the list contains real data)
#'
#' @export
#'

perm_mat <- function(X=NULL, perms=NULL, columns=FALSE){
	
	X_list <- vector("list", length = length(perms))
	X_list <- lapply(X_list, function(x) X)
	
	for(i in 2:length(X_list)){ #assign the new rownames, the first element are real data
		rownames(X_list[[i]]) <- perms[[i]]
		if(columns){
			colnames(X_list[[i]]) <- perms[[i]]
		}
	}
	
	#reorder the rownames (and optinally the columns) of the permuted X_i as in the original matrix (the first element)
	X_list <- lapply(X_list, function(X_i) X_i[match(rownames(X_list[[1]]), rownames(X_i)), , drop = F])
	if(columns){
		X_list <- lapply(X_list, function(X_i) X_i[, match(colnames(X_list[[1]]), colnames(X_i)), drop = F])
	}
	
	return(X_list)
	
}
