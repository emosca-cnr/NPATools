#' Assess modularity
#' @param G igraph object
#' @param vertices set of (top ranking) vertices to study
#' @param ranks ranks at which the modularity will be computes
#' @param commMethods community method(s) that will be considered
#' @param BPPARAM BBPARAM
#' @param minComponentSize minimum size of connected components
#' @export
#' @import BiocParallel igraph
#' 
assess_modularity <- function(G = NULL, vertices = NULL, ranks = NULL, commMethods=c("fastgreedy", "multilev"),  BPPARAM = NULL, minComponentSize=2) {
  
  if (is.null(ranks)) {
    ranks <- seq(10, length(vertices), by = 10)
  }
  
  if (is.null(BPPARAM)) {
    BPPARAM <- SerialParam()
  }
  
  top_nets <- lapply(ranks, function(k) get_nconn_comp(induced_subgraph(graph = G, vids = V(G)$name[V(G)$name %in% vertices[1:k]]), n = minComponentSize))
  names(top_nets) <- ranks
  
  #top_nets_conn <- lapply(ranks, function(k) get_max_conn_comp(induced_subgraph(graph = G, vids = V(G)$name[V(G)$name %in% vertices[1:k]])))
  top_nets_conn <- lapply(ranks, function(k) largest_component(induced_subgraph(graph = G, vids = V(G)$name[V(G)$name %in% vertices[1:k]])))
  names(top_nets_conn) <- ranks
  
  
  top_comm <- BiocParallel::bplapply(top_nets, function(g_i) {
    #print(x)
    comm <- find_communities(g_i, verbose = F, methods = commMethods)
    best <- comm$info
    best <- best[order(-best$modularity, best$n),]
    tmp <- array(c(vcount(g_i), ecount(g_i), best[1, ]))
    return(tmp)
  }, BPPARAM = BPPARAM)
  top_comm <- do.call(rbind, top_comm)
  
  top_comm_conn <- BiocParallel::bplapply(top_nets_conn, function(g_i) {
    #print(x)
    comm <- find_communities(g_i, verbose = F, methods = commMethods)
    best <- comm$info
    best <- best[order(-best$modularity, best$n),]
    tmp <- array(c(vcount(g_i), ecount(g_i), best[1, ]))
    return(tmp)
  }, BPPARAM = BPPARAM)
  top_comm_conn <- do.call(rbind, top_comm_conn)
  
  df.out <- data.frame(rank = names(top_nets),
                       n_vertex = as.numeric(top_comm[, 1]),
                       n_edge = as.numeric(top_comm[, 2]),
                       algorithm = unlist(top_comm[, 3]),
                       modularity = as.numeric(top_comm[, 4]),
                       n_community = as.numeric(top_comm[, 5]),
                       n_vertex_max = as.numeric(top_comm_conn[, 1]),
                       n_edge_max = as.numeric(top_comm_conn[, 2]),
                       algorithm_max = unlist(top_comm_conn[, 3]),
                       modularity_max = as.numeric(top_comm_conn[, 4]),
                       n_community_max = as.numeric(top_comm_conn[, 5]),
                       stringsAsFactors = F)
  
  return(df.out)
}