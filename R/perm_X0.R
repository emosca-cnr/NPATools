#' Permutation of X0
#'
#' Permutation of input matrix X0
#' @param X0 matrix; a score matrix X0, in which each column (layer) is a score vector over all vertices of A.
#' @param perms list of permutations of the rownames of X0, obtained by means of `permute_vertices()`
#' @return list with permutations of X0 (where the first element of the list contains real data)
#'
#' @export
#' @importFrom Matrix Matrix
#' @importFrom ggplot2 cut_interval cut_number cut_width
#' @importFrom stats setNames
#'

perm_X0 <- function(X0=NULL, perms=NULL){
	
	X0_list <- vector("list", length = length(perms))
	X0_list <- lapply(X0_list, function(x) X0)
	
	for(i in 2:length(X0_list)){
		rownames(X0_list[[i]]) <- perms[[i]]
	}
	
	X0_list <- lapply(X0_list, function(X0_i) X0_i[match(rownames(X0_list[[1]]), rownames(X0_i)), , drop = F])
	
	return(X0_list)
	
}
