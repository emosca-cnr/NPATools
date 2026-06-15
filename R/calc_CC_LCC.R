#' Calculate connected components and largest connected componets
#' @param G igraph object
#' @param vertices set of (top ranking) vertices to study
#' @param ranks optional ranks at which the modularity will be computes
#' @param vertices_max maximum number of vertices allowed; the higher vertices_max the higher the computational cost
#' @export
#' @importFrom igraph induced_subgraph largest_component


calc_CC_LCC <- function(G = NULL, vertices = NULL, ranks=NULL, vertices_max=500) {
	
	stopifnot(!is.null(G) | !is.null(vertices))
	
	stopifnot(length(vertices) <= vertices_max)
	
	if (is.null(ranks)) {
		ranks <- seq(5, length(vertices), by = 5)
	}
	
	###CC and LCC
	g_CC <- lapply(ranks, function(k) get_nconn_comp(induced_subgraph(graph = G, vids = V(G)$name[V(G)$name %in% vertices[1:k]]), n = 2))
	names(g_CC) <- ranks
	
	g_LCC <- lapply(ranks, function(k) largest_component(induced_subgraph(graph = G, vids = V(G)$name[V(G)$name %in% vertices[1:k]])))
	names(g_LCC) <- ranks

	return(list(g_CC=g_CC, g_LCC=g_LCC))
	
}

