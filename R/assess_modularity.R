#' Assess modularity
#' @param g_CC_LCC list of CC and LCC as produced by means of calc_CC_LCC()
#' @param commMethods community method(s) that will be considered
#' @param BPPARAM BBPARAM
#' @export
#' @importFrom  BiocParallel bplapply SerialParam
#' @importFrom igraph vcount ecount

assess_modularity <- function(g_CC_LCC = NULL, commMethods=c("fastgreedy", "multilev"), BPPARAM = NULL) {
  
  if (is.null(BPPARAM)) {
    BPPARAM <- SerialParam()
  }
  
  top_comm <- bplapply(g_CC_LCC$g_CC, function(g_i) {
    #print(x)
    comm <- find_communities(g_i, verbose = F, methods = commMethods)
    best <- comm$info
    best <- best[order(-best$modularity, best$n),]
    tmp <- array(c(vcount(g_i), ecount(g_i), best[1, ]))
    return(tmp)
  }, BPPARAM = BPPARAM)
  
  top_comm <- do.call(rbind, top_comm)
  
  top_comm_conn <- bplapply(g_CC_LCC$g_LCC, function(g_i) {
    #print(x)
    comm <- find_communities(g_i, verbose = F, methods = commMethods)
    best <- comm$info
    best <- best[order(-best$modularity, best$n),]
    tmp <- array(c(vcount(g_i), ecount(g_i), best[1, ]))
    return(tmp)
  }, BPPARAM = BPPARAM)
  top_comm_conn <- do.call(rbind, top_comm_conn)
  
  df.out <- data.frame(rank = names(g_CC_LCC$g_CC),
                       n_V_CC = as.numeric(top_comm[, 1]),
                       n_E_CC = as.numeric(top_comm[, 2]),
                       a_CC = unlist(top_comm[, 3]),
                       Q_CC = as.numeric(top_comm[, 4]),
                       n_comm_CC = as.numeric(top_comm[, 5]),
                       n_V_LCC = as.numeric(top_comm_conn[, 1]),
                       n_E_LCC = as.numeric(top_comm_conn[, 2]),
                       a_LCC = unlist(top_comm_conn[, 3]),
                       Q_LCC = as.numeric(top_comm_conn[, 4]),
                       n_comm_LCC = as.numeric(top_comm_conn[, 5]),
                       stringsAsFactors = F)
  
  return(df.out)
}