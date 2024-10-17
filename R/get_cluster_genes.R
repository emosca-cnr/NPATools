#' Get pathway clusters from enrichment map
#' @param gene_set_clusters gene set cluster data frame created by get_gene_set_clusters()
#' @param gene_set_list gene set list
#' @export

get_cluster_genes <- function(gene_set_clusters=NULL, gene_set_list=NULL){
  
  return(tapply(gene_set_clusters$name, gene_set_clusters$comm_id, function(x) sort(as.character(unique(unlist(gene_set_list[x]))))))
  
  
}