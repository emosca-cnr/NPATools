#' pc_wmz Calculate participation coefficient and within module degree z-score
#'
#' @description Calculate participation coefficient
#' and within module degree z-score
#' @param topNetwork network with some parameters: "comm_id" for community
#' membership and "label" for labels (for plotting purpose).
#' @param labelBy vertex attribute of topNetwork to use for labels
#' @param useP whether to use p or not in the vertical axis
#' @param hCut value to divide the horizontal axis
#' @return dataframe with pc, wmz, region, description and role
#' @importFrom brainGraph part_coeff
#' @importFrom igraph V get.vertex.attribute
#' @export

functional_cartography <- function(topNetwork = NULL, labelBy = NULL, useP = FALSE, hCut = NULL) {
  
  stopifnot(length(V(topNetwork)$comm_id) > 0)
  
  # partecipation coefficient
  pc <- part_coeff(topNetwork, V(topNetwork)$comm_id)
  
  # within module z-score
  #wmZ <-brainGraph::within_module_deg_z_score(topNetwork,
  #                                              V(topNetwork)$comm_id)
  
  wmZ <-
    within_module_degree(g = topNetwork,
                         memb = as.numeric(V(topNetwork)$comm_id),
                         useP = useP)
  #wm_z[is.infinite(wm_z[, 1]), 1] <- 0
  wmZ[is.na(wmZ[, 1]), 1] <- 0
  wmZ[is.na(wmZ)] <- 0
  
  df <-
    data.frame(
      name = V(topNetwork)$name,
      comm_id = V(topNetwork)$comm_id,
      P = as.numeric(pc),
      wmZ,
      R = NA,
      definition = NA,
      stringsAsFactors = F
    )
  
  wmZ <- wmZ[, 1]
  maxZscore <- max(wmZ)
  minZscore <- min(wmZ - 1)
  
  if (is.null(hCut)) {
    if (useP) {
      hCut <- -log10(0.05)
    } else {
      hCut <- (maxZscore * 2.5) / 8
    }
  }
  
  # look for coordinates
  # NON-HUBS : z < 2.5
  # for R1: ULTRA PERIPHERAL NODES (p<=0.05) = nodes with all
  #their links within their module
  df$R[df$P <= 0.05 & df$wmd_score < hCut] <- "R1"
  df$definition[df$R == "R1"] <- "ULTRA PERIPHERAL NODES"
  
  # for R2: PERIPHERAL NODES (0.05 < p <= 0.62) = nodes with
  # most links within their module
  df$R[df$P > 0.05 &
         df$P <= 0.62 & df$wmd_score < hCut] <- "R2"
  df$definition[df$R == "R2"] <- "PERIPHERAL NODES"
  
  # for R3: NON-HUB CONNECTOR NODES (0.62 < p <= 0.8) = nodes
  # with many links to other modules
  df$R[df$P > 0.62 &
         df$P <= 0.8 & df$wmd_score < hCut] <- "R3"
  df$definition[df$R == "R3"] <- "NON-HUB CONNECTOR NODES"
  
  # for R4: NON-HUB KINLESS NODES (p > 0.8) = nodes with links
  # homogeneously distributed among all modules
  df$R[df$P > 0.8 & df$wmd_score < hCut] <- "R4"
  df$definition[df$R == "R4"] <- "NON-HUB KINLESS NODES"
  
  ### HUBS : z > 2.5
  # for R5: PROVINCIAL HUBS (p <= 0.3) = hub nodes with the vast
  # majority of links within their module
  df$R[df$P <= 0.3 & df$wmd_score > hCut] <- "R5"
  df$definition[df$R == "R5"] <- "PROVINCIAL HUBS"
  
  # for R6: CONNECTOR HUBS (0.3 < p <= 0.75) = hubs with many links
  # to most of the other modules
  df$R[df$P > 0.3 &
         df$P <= 0.75 & df$wmd_score > hCut] <- "R6"
  df$definition[df$R == "R6"] <- "CONNECTOR HUBS"
  
  # for R7: KINLESS HUBS (p > 0.75) = hubs with links homogeneously
  # distributed among all modules
  df$R[df$P > 0.75 & df$wmd_score > hCut] <- "R7"
  df$definition[df$R == "R7"] <- "KINLESS HUBS"
  
  if (!is.null(labelBy)) {
    vlab <- get.vertex.attribute(topNetwork, name = labelBy)
    df <- cbind(df, vlab)
  }
  
  V(topNetwork)$R <- df$R
  V(topNetwork)$P <- df$P
  V(topNetwork)$wmd_score <- df$wmd_score
  
  
  return(topNetwork)
  
}
