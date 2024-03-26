#' Assign communities
#' @description the community will be added as the vertex attribute "comm_id"
#' @param g igraph object
#' @param commList community structure
#' @export
#' 

assign_communities <- function(g=NULL, commList=NULL){
  
  V(g)$comm_id <- commList$membership[match(V(g)$name, commList$names)]
  #V(g)$comm_id <- commList$membership
  print(table(V(g)$comm_id))
  
  return(g)
  
}