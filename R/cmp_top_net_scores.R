#' Compare top networks scores
#'
#' @description Define significance for enriched genes
#' @param netEnrRes dataframe resulting from `assess_enrichment()`
#' @param NRRes dataframe resulting from network resampling from package dmfind `dmfind::NR()`
#' @param top number of top ranks that will be plotted in red
#' @param out_dir output directory
#' @export
#' @importFrom plotrix thigmophobe.labels
#' @importFrom graphics par points
#' @importFrom grDevices png dev.off

cmp_top_net_scores <-
  function(NRRes = NULL,
           netEnrRes = NULL,
           top=5,
           out_dir=NULL) {
    
    if(is.null(out_dir)){
      out_dir <- getwd()
    }else{
      dir.create(out_dir, recursive = T, showWarnings = F)
    }

    ans <-
      merge(
        NRRes$NRsummary[, c("rank", "p_val")],
        netEnrRes$en_summary[, c("rank", "p_val")],
        by = "rank",
        suffixes = c("_NR", "_NE"),
        sort = F
      )
    ans <- ans[order(as.numeric(ans$rank)),]
    
    xx <- ans[, 2]
    yy <- ans[, 3]
    
    #xyarea <- xx_norm*yy_norm
    #ans$score <- xx + yy
    ans$score <- xx * yy
    ans$score <- ans$score - (ans$score * log(ans$score))
    ans$score <- -log10(ans$score)
    
    xx <- -log10(xx)
    yy <- -log10(yy)
    
    idx_max <- order(-ans$score)[1:top]
    cat("Smallest highest scoring network at rank: ", ans$rank[idx_max[1]], "\n")
        
    png(file.path(out_dir, "NR_vs_NE.png"), width = 180, height = 100, units="mm", res=300)
    
    par(mfrow = c(1, 2))
    par(mar = c(3, 3, 1, 1))
    par(mgp = c(1.5, 0.5, 0))
    
    
    if (nrow(unique(ans[, 2:3])) != nrow(ans)) {
      cat("Overlapping points detected, adding some jitter\n")
      xx <- jitter(xx)
      yy <- jitter(yy)
    }
    
    plot(
      xx,
      yy,
      xlab = paste0("p1 (NR)"),
      ylab = paste0("p2 (NE)"),
      pch = 16,
      cex = 0.6,
      cex.axis = 0.7,
      cex.lab = 0.7
    )
    points(xx[idx_max],
           yy[idx_max],
           pch = 16,
           col = "red",
           cex = 0.6)
    
    text_col <- rep("black", length(xx))
    text_col[idx_max] <- "red"
    thigmophobe.labels(xx, yy, ans$rank, cex = 0.6, col = text_col)
    
    plot(
      as.numeric(ans$rank),
      ans$score,
      xlab = "rank",
      ylab = "-log10( P(P1*P2 < p1*p2) )",
      pch = 16,
      cex = 0.6,
      cex.axis = 0.7,
      cex.lab = 0.7,
      col=text_col
    )
    dev.off()
    
    return(ans)
    
  }
