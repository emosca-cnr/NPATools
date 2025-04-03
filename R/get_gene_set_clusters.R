#' Get pathway clusters from enrichment map
#' @param enrichment_map_graph igraph object created by enrichment_map()
#' @param gene_set_list gene set list
#' @export
#' @importFrom igraph vertex_attr


get_gene_set_clusters <- function(enrichment_map_graph=NULL, gene_set_list=NULL){
  
  ans <- as.data.frame(vertex_attr(enrichment_map_graph))
  ans$genes <- sapply(ans$name, function(x) paste0(sort(gene_set_list[[x]]), collapse = ";"))
  
  return(ans)
  
}