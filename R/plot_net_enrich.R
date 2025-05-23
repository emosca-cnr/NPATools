#' Plot network enrichment results
#'
#' @param netEnrRes output of `assess_enrichment()`
#' @param showTopSig if not NULL, an integer that indicatesd the number of top significant ranks to display
#' @param sigStat the quantity to be plotted: for ORA, can be any of "er", "p" or "p_adj"; for GSEA can be "es", "p", "p_val", "nes" or "FDRq"
#' @export
#' @importFrom graphics lines

plot_net_enrich <-
  function(netEnrRes = NULL,
           showTopSig = 10,
           sigStat = "p_val"){
    
    en_res <- netEnrRes$en_summary[order(as.numeric(netEnrRes$en_summary$rank)), ]

    par(mgp = c(1.5, 0.5, 0))
    par(mar = c(3, 3, 3, 1))
    
    if (sigStat %in% c("p_val", "p_adj", "FDRq")) {
      yy <- log10(en_res[, sigStat])
      ylab <- paste0("log10(", sigStat, ")")
    } else{
      yy <- en_res[, sigStat]
      ylab <- sigStat
    }
    
    top_10_idx <-
      order(yy, decreasing = !(sigStat %in% c("p_val", "p_adj", "FDRq")))[1:showTopSig]
    
    plot(
      en_res$rank,
      yy,
      pch = 16,
      xlab = "rank",
      ylab = ylab,
      cex = 0.6,
      main = "enrichment"
    )
    lines(en_res$rank, yy)
    
    if (!is.null(showTopSig)) {
      points(en_res$rank[top_10_idx],
             yy[top_10_idx],
             pch = 16,
             cex = 0.6,
             col = "red")
      thigmophobe.labels(en_res$rank[top_10_idx], yy[top_10_idx], en_res$rank[top_10_idx], cex =
                           0.7)#, pos = 3)
    }
    
}
