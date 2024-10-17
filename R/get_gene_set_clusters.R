#' Get pathway clusters from enrichment map
#' @param enrichment_map_graph igraph object created by enrichment_map()
#' @export
#' @importFrom igraph vertex_attr


get_gene_set_clusters <- function(enrichment_map_graph=NULL){
  
  return(as.data.frame(vertex_attr(enrichment_map_graph)))
  
}