#' Get pathway clusters from enrichment map
#' @param gene_set_clusters gene set cluster data frame created by get_gene_set_clusters()
#' @param gene_set_list gene set list
#' @export

get_cluster_genes <- function(gene_set_clusters=NULL, gene_set_list=NULL){
  
  clust_genes <- tapply(gene_set_clusters$name, gene_set_clusters$comm_id, function(x) sort(as.character(unique(unlist(gene_set_list[x])))))
  
  clust_genes_df <- data.frame(clust_id=rep(names(clust_genes), lengths(clust_genes)), gene_id=unlist(clust_genes), size=rep(lengths(clust_genes), lengths(clust_genes)))
  clust_genes_df$cg_id <-apply(clust_genes_df[, 1:2], 1, paste0, collapse="_")

  return(clust_genes_df)
  
}