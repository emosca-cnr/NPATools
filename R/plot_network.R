#' Plot a network
#' 
#' @description Extraction and plotting of connected subnetworks 
#' from a graph and a set of selected vertices
#' @param graph igraph object
#' @param colorBy object to be colored by an attribute
#' @param colorQuant TRUE/FALSE
#' @param labelBy select a feature to label the graph
#' @param pal palette to color the graph
#' @param community communities from network analysis 
#' @param commWin distance of the vertices inside a community
#' @param commWb distance between communities
#' @param lo layout 
#' @param legendOff legend
#' @param ... additional arguments to `plot.igraph()`
#' @export
#' @import igraph
#' @importFrom plotrix thigmophobe.labels

plot_network <- function(graph=NULL, colorBy=NULL, colorQuant=TRUE, 
                         labelBy="name", pal=NULL, community=NULL, commWin=2, commWb=1, lo=NULL, 
                         legendOff=FALSE, labelComm=TRUE, labelCommCex=0.8, labelCommCol="black", ...){
  
  if(!is.null(colorBy)){
    colorValues <- get.vertex.attribute(graph, colorBy)
    if(colorQuant){
      colorValues <- cut(colorValues, length(pal), dig.lab = 1)
    }
    V(graph)$color <- pal[as.numeric(as.factor(colorValues))]
  }
  if(labelBy == ""){
    vertexLabels <- ""
  }else{
    vertexLabels <- get.vertex.attribute(graph, labelBy)
  }
  V(graph)$label <- vertexLabels
  
  if(labelComm & !length(V(graph)$comm_id)>0){
    cat("To show community labels, please define V(graph)$comm_id.\n")
    labelComm <- FALSE
  }
  
  ew <- NULL
  if(!is.null(community)){
    ew <- edge_weights(community, graph, weightWithin=commWin,
                       weightBetween=commWb)
  }
  
  return_lo <- FALSE
  if(is.null(lo)){
    lo <- layout_with_fr(graph, weights=ew)
    return_lo <- TRUE
  }

  if(labelComm){
    network_ann_coords <- merge(tapply(lo[, 1], INDEX = V(graph)$comm_id, mean), tapply(lo[, 2], INDEX = V(graph)$comm_id, mean), by=0) 
    comm_lab_xy <- norm_coords(as.matrix(network_ann_coords[, 2:3]), -1, 1, -1, 1)
  }
  
    
  do_lgd <- FALSE
  if(!is.null(colorBy) & !legendOff){
    if(colorBy != "comm_id"){
      do_lgd <- TRUE
      layout(matrix(c(1, 2), ncol=2), widths = c(0.9, 0.2))
    }
  }
  par(mar=c(1, 1, 1, 1))
  
  plot.igraph(graph, vertex.labels=vertexLabels, layout=lo, ...)
  if(labelComm){
    thigmophobe.labels(comm_lab_xy[, 1], comm_lab_xy[, 2], labels = network_ann_coords[, 1], cex=labelCommCex, font=2, col=labelCommCol)
  }
  
  if(do_lgd){
    plot.new()
    legend("bottomright", levels(factor(colorValues, 
                                        levels=sort(unique(colorValues)))), 
           col=pal, pch=16, bty="n", cex=0.8)
  }
  
  if(return_lo){
    return(lo)
  }
  
}
