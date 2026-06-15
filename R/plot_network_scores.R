#' Plot the network scoring summary
#'
#' @description Plot the values obtained by means of score_networks_summary()
#' @param net_scores_summm the table provided by the function score_networks_summary()
#' @importFrom ggplot2 ggplot aes geom_line set_theme scale_color_manual geom_col  geom_point geom_segment labs geom_hline
#' @importFrom reshape2 melt
#' @importFrom pals brewer.paired
#' @returns A ggplot2 object
#' @export

plot_network_scores <- function(net_scores_summm = NULL){
  
  rank <- value <- variable <- peak <- NULL # Setting the variables to NULL first
  
  set_theme(theme_science())
  col <- brewer.paired(5)
  
  net_scores_summm$is.max <- FALSE
  net_scores_summm$is.max[which.max(net_scores_summm$y)] <- TRUE
  
  plot_data <- melt(net_scores_summm[, 1:6], id.vars = "rank")
  plot_data$ptype <- NA
  plot_data$ptype[plot_data$rank %in% net_scores_summm$rank[net_scores_summm$peak]] <- 16
  plot_data$vltype <- NA
  plot_data$vltype[plot_data$rank %in% net_scores_summm$rank[net_scores_summm$is.max]] <- 2
  
  p <- ggplot(plot_data, aes(x = as.numeric(rank), y = value, color = variable)) +
    geom_line() +
    geom_point(data=plot_data[plot_data$variable == "y", ], shape = plot_data$ptype[which(plot_data$variable == "y")]) +
    geom_segment(data=plot_data[plot_data$variable == "y", ], aes(y=0, yend=plot_data$value[plot_data$variable == "y"]), lty=plot_data$vltype[which(plot_data$variable == "y")]) +
    geom_hline(yintercept = 0, lty=2) +
    scale_color_manual(values = col) +
    labs(x = "rank", y="", color="")
  
  return(p)
  
}
