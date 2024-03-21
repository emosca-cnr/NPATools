#' Enrichment map
#' @description This functions calculates similarities between gene-sets and plot a resulting enrichment map
#' @details enrichment_map() function calculates similarities between each gene-set pair by using `method` metric. Subsequently,
#'  these similarities are filtered to maintain the ones >= `coeff`. `comm_method` algorithm is then used to identify communities, 
#'  which may be filtered to plot only the ones composed by at least `min_comm_size` gene-sets. Enrichment map is then plotted by using
#'  `ggraph` package.
#' @param gs_score named vector of pathway scores
#' @param coeff threshold for the similarity score between two gene sets
#' @param min_comm_size minimum size of a community to be considered in the enrichment map. If `min_comm_size = 1` and `all_gs = TRUE`,
#'  then all gene-sets are displayed
#' @param gs_size named vector with size f the gene-sets in `gs_list`
#' @param set_sim_df optional, data.frame with three columns 'set1', 'set2' and 'sim'
#' @return The function returns a list containing: 
#' \itemize{
#'  \item igraph = network object used for plotting the enrichment map
#'  \item network_data = data.frame with layout (X1,X2 columns), name of the gene-sets together with community ids ("comm_id"), 
#'   score and number of genes ("n_genes")
#'  \item path_comm_genes = list composed by the genes present in each community
#'  \item sim_coeff = data.frame with the similarity score calculated between each gene-set pair
#'  \item plot = enrichment map plot obtained by using `ggraph` package functions. Only if `save_plot = FALSE`
#' }
#' @import igraph
#' @export
#'
enrichment_map <- function(gs_scores=NULL, gene_set_sim=NULL, min_sim=0.7, clustering_f=cluster_fast_greedy, min_comp_size=1,
                           gs_size=NULL, add.top.k=NULL){
  #method=c('overlap')
  #similarity method and coeff----------------
  
  gene_set_sim <- as.matrix(gene_set_sim)
  gene_set_sim[gene_set_sim < min_sim] <- 0
  
  #building pathway network + communities------------------------------
  path_g <- graph_from_adjacency_matrix(adjmatrix = gene_set_sim, mode = "undirected", weighted = TRUE, diag = FALSE) 
  
  #remove components
  if(min_comp_size>1){
    
    comp <- clusters(path_g)
    clstrOk <- which(comp$csize >= min_comp_size)
    nodesOk <- which(comp$membership %in% clstrOk)
    path_g <- induced.subgraph(path_g, nodesOk)
  }
  
  print(path_g)
  
  cat("Community detection...\n")
  path_g_clusters <- clustering_f(path_g, weights = E(path_g)$weight)
  cat("modularity =", modularity(path_g_clusters), "\n")
  cat("cluster sizes =", sizes(path_g_clusters), "\n")

  V(path_g)$comm_id <- as.numeric(membership(path_g_clusters))
  
  V(path_g)$score <- gs_scores[match(V(path_g)$name, names(gs_scores))] #pathway score based on ORA res
  V(path_g)$gs_size <- gs_size[match(V(path_g)$name, names(gs_size))]
  
  return(path_g)
 
}

