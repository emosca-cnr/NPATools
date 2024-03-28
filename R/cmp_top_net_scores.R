#' Compare top networks scores
#'
#' @description Define significance for enriched genes
#' @param netEnrRes dataframe resulting from `assess_enrichment()`
#' @param NRRes dataframe resulting from network resampling from package dmfind `dmfind::NR()`
#' @param columm which column to consider; in case multiple runs have been done
#' @param sigStatNR which p-value will be considered from NRRes
#' @param sigStatEn which p-value will be considered from netEnrRes
#' @param norm if TRUE (default) the two disitribution of -log10(p) will be mapped in the interval [1, 2]
#' @param top number of top ranks that will be plotted in red
#' @export
#' @importFrom plotrix thigmophobe.labels
#'

cmp_top_net_scores <-
  function(NRRes = NULL,
           netEnrRes = NULL,
           column = 1,
           sigStatNR = "p",
           sigStatEn = "p",
           norm = FALSE,
           top=5) {

    ans <-
      merge(
        NRRes$NRsummary[, c("rank", sigStatNR)],
        netEnrRes$en_summary[, c("rank", sigStatEn)],
        by = "rank",
        suffixes = c("_NR", "_NE"),
        sort = F
      )
    ans <- ans[order(as.numeric(ans$rank)),]
    
    xx <- -log10(ans[, 2])
    #yy <- sapply(as.numeric(netEnrRes$en_summary$id), function(x) sum(p_sorted[1:x] < sigThr) / x)
    yy <- -log10(ans[, 3])
    
    if (norm) {
      xx <- 1 + (xx - min(xx)) / (max(xx) - min(xx))
      yy <- 1 + (yy - min(yy)) / (max(yy) - min(yy))
    }
    
    #xyarea <- xx_norm*yy_norm
    ans$score <- xx * yy
    
    par(mfrow = c(1, 2))
    par(mar = c(3, 3, 1, 1))
    par(mgp = c(1.5, 0.5, 0))
    
    idx_max <- order(-ans$score)[1:top]
    
    if (nrow(unique(ans[, 2:3])) != nrow(ans)) {
      cat("Overlapping points detected, adding some jitter\n")
      xx <- jitter(xx)
      yy <- jitter(yy)
    }
    
    plot(
      xx,
      yy,
      xlab = paste0("NR [-log10(", sigStatNR, ")]"),
      ylab = paste0("NE [-log10(", sigStatEn, ")]"),
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
      ylab = paste0("NR[-log10(", sigStatNR, ")] * NE[-log10(", sigStatEn, ")]"),
      pch = 16,
      cex = 0.6,
      cex.axis = 0.7,
      cex.lab = 0.7,
      col=text_col
    )
    
    return(ans)
    
  }
