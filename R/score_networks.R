#' Score the networks composed of a ranked list of genes
#' 
#' @description
#' The networks are scored by:
#' 
#' * modularity;
#' * network components;
#' * omega z-score peaks;
#' * presence of linkers.
#' @returns
#' A data frame with the following quantities for CC and LCC:
#' 
#' * `X0_`: sum of the X0 values in the network divided by the sum of the n genes with the highest X0 values, where n network size (CC or LCC);
#' * `linkers`: fraction of linkers, that is, genes with X0(i) equal to 0;
#' * `omega`: NR omega value;
#' * `z`: omega z-score;
#' * `p_val`: omega p_val;
#' * `n_V_`: number of vertices;
#' * `n_E_`: number of links;
#' * `a_`: algorithm used to find the communities;
#' * `Q_`: modularity;
#' * `n_comm_`: number of communities.
#' @param g_CC_LCC list of CC and LCC as produced by means of calc_CC_LCC()
#' @param NRsummary NRsummary element of the output of function NR()
#' @param X0 the X0 matrix
#' @param X0_column the column of X0 that will be considered
#' @param commMethods community method(s) that will be considered
#' @param BPPARAM optional BBPARAM
#' @export
#' @importFrom igraph V vcount
#' @importFrom  BiocParallel SerialParam
#' @importFrom stats setNames

score_networks <- function(g_CC_LCC = NULL, X0=NULL, X0_column=1, NRsummary=NULL, commMethods="fastgreedy", BPPARAM=NULL) {
	
	stopifnot(!is.null(g_CC_LCC) | !is.null(NRsummary) | !is.null(X0))
	
	if (is.null(BPPARAM)) {
		BPPARAM <- SerialParam()
	}
	

	cat("Quantyfing the X0 enrichment...\n")
	cat("\tX0_column =", X0_column, ",", colnames(X0)[1], "...\n")
	X0_vector <- setNames(X0[, X0_column], rownames(X0))
	X0_vector_cumsum <- cumsum(sort(X0_vector, decreasing = T)) #cumulative sum of of gene scores in decreasing order
	
	X0_sum_CC <- unlist(lapply(g_CC_LCC$g_CC, function(x) as.numeric(sum(X0_vector[names(X0_vector) %in% V(x)$name]) / X0_vector_cumsum[vcount(x)]))) #scores in the CC relative the maximum score that can be achieved by a CC of the same size
	
	X0_sum_LCC <- unlist(lapply(g_CC_LCC$g_LCC, function(x) as.numeric(sum(X0_vector[names(X0_vector) %in% V(x)$name]) / X0_vector_cumsum[vcount(x)])))#scores in the LCC relative the maximum score that can be achieved by a CC of the same size
	
	cat("Qauntifying the presence of linkers...\n")
	linkers <- unlist(lapply(g_CC_LCC$g_CC, function(x) as.numeric(sum(X0_vector[names(X0_vector) %in% V(x)$name]==0) / vcount(x))))

	cat("Assessing modularity...\n")
	mod_res <- assess_modularity(g_CC_LCC = g_CC_LCC, commMethods = commMethods, BPPARAM = BPPARAM)
	
	cat("Assembling output table...\n")
	ans <- data.frame(rank=names(X0_sum_CC), X0_CC=X0_sum_CC, X0_LCC=X0_sum_LCC, linkers=linkers)
	
	ans <- merge(ans, NRsummary[, c("rank", "omega", "z", "p_val")], by="rank", sort=F)
	ans <- merge(ans, mod_res[, !grepl("algorithm", colnames(mod_res))], by="rank", sort=F)

	ans <- ans[order(as.numeric(ans$rank)), ]
	
	return(ans)
	
}

