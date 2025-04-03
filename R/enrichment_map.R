#' Enrichment map
#' @description This functions calculates similarities between gene-sets and plot a resulting enrichment map
#' @details enrichment_map() function calculates similarities between each gene-set pair by using `method` metric. Subsequently,
#'  these similarities are filtered to maintain the ones >= `min_sim`. `clustering_f` function is then used to identify communities, 
#'  which may be filtered to plot only the ones composed by at least `min_comp_size` gene-sets. 
#' @param gs_scores named vector of pathway scores
#' @param gene_set_sim gene set similarity calculated through `calc_gs_sim()`
#' @param min_sim threshold for the similarity score between two gene sets
#' @param clustering_f graph clustering function, `igraph::cluster_fast_greedy()` by default 
#' @param min_comp_size minimum size of a community to be considered in the enrichment map
#'  then all gene-sets are displayed
#' @param gs_size named vector with size of the gene-sets in `gs_scores`
#' @param verbose verbosity
#' @importFrom igraph graph_from_adjacency_matrix clusters induced.subgraph V cluster_fast_greedy E V<- membership sizes
#' @export
#'
enrichment_map <- function(gs_scores=NULL, gene_set_sim=NULL, min_sim=0.2, clustering_f=cluster_fast_greedy, min_comp_size=1, gs_size=NULL, verbose=TRUE){
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
  
  if(verbose){
    cat("Community detection...\n")
  }
  path_g_clusters <- clustering_f(path_g, weights = E(path_g)$weight)
  
  Q <- modularity(path_g_clusters)
  if(verbose){
    cat("modularity =", Q, "\n")
  }
  
  cl_sizes <- sizes(path_g_clusters)

  if(verbose){
    cat("Mean cluster size =", mean(cl_sizes), "\n")
  }
  if(verbose){
    cat("# clusters size <= 2", sum(cl_sizes[cl_sizes<=2]), "\n")
  }
  V(path_g)$comm_id <- as.numeric(membership(path_g_clusters))
  
  V(path_g)$score <- gs_scores[match(V(path_g)$name, names(gs_scores))] #pathway score based on ORA res
  V(path_g)$gs_size <- gs_size[match(V(path_g)$name, names(gs_size))]
  
  return(list(ig=path_g, Q=Q, N=length(cl_sizes), avg_size=mean(cl_sizes), comp_lt2=sum(cl_sizes[cl_sizes<=2])))
 
}

