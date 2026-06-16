#' Plot the functional cartography of the network
#' @param topNetworks igraph object
#' @param labelBy optional, one of the vertex attributes
#' @param hCut optional, where to place the division of the vertical axis
#' @param useP whether the p-value associated with the z-score was used
#' @export
#' @importFrom igraph vcount vertex_attr
#' @importFrom ggplot2 set_theme ggplot coord_cartesian geom_rect scale_fill_manual geom_point geom_text theme
#' @returns ggplo2 object
 
plot_fc <- function(topNetworks = NULL, labelBy = NULL, hCut = NULL, useP = FALSE){
  
  x <- xend <- y <- yend <- class <- label <- P <- z <- NULL # Setting the variables to NULL first
  
  set_theme(theme_science())
  
  if (!is.null(labelBy)) {
    vlab <- vertex_attr(topNetworks, name = labelBy)
  }else{
    vlab <- rep("", vcount(topNetworks))
  }
  
  plot_data <- data.frame(P=V(topNetworks)$P, z=V(topNetworks)$wmd_score, label=vlab)
  
  maxZscore <- max(plot_data$z)
  minZscore <- min(plot_data$z)
  
  if (is.null(hCut)) {
    if (useP) {
      hCut <- -log10(0.05)
    } else {
      hCut <- (maxZscore * 2.5) / 8
    }
  }
  
  plot_rects <- data.frame(x=c(-1, 0.05, 0.62, 0.8, -1, 0.3, 0.75), xend=c(0.05, 0.62, 0.8, 2, 0.3, 0.75, 2), y=c(minZscore-1, minZscore-1, minZscore-1, minZscore-1, hCut, hCut, hCut), yend=c(hCut, hCut, hCut, hCut, maxZscore+1, maxZscore+1, maxZscore+1), class=c("R1", "R2", "R3", "R4", "R5", "R6", "R7"))
  
  col <- c("ivory4", "lightsalmon", "darkseagreen", "slateblue2", "lightyellow", "mistyrose", "gray90")
  
  p <- ggplot() +
    coord_cartesian(xlim = c(0, 1), ylim = c(minZscore, maxZscore)) +
    geom_rect(data = plot_rects, aes(xmin = x, xmax = xend, ymin = y, ymax = yend, fill = class)) +
    scale_fill_manual(values = col) +
    geom_point(data = plot_data, aes(x=P, y=z)) +
    geom_text(data = plot_data, aes(x=P, y=z, label=label), nudge_y = -0.1) +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
  
  return(p)
  
}

