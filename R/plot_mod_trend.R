#' Plot modularity trend
#' @param assessModRes result of `assess_modularity()`
#' @param file filename
#' @return plot
#' @examples
#' \dontrun{plot_S(nsiRes=NULL, X0, initSigGenes=NULL, file="S.jpg")}
#' @usage plot_S(nsiRes=NULL, X0, initSigGenes=NULL, file="S.jpg")
#' @export

plot_modu_trend <- function(assessModRes = NULL) {
        
        layout(matrix(1:2, ncol = 2), widths = c(0.9, 0.1))
        
        par(mar = c(2.5, 2.5, 1, 2))
        par(mgp = (c(1.5, .5, 0)))
        
        plot(
            assessModRes$rank,
            assessModRes$modularity,
            xlab = "rank",
            ylab = "Q",
            pch = 16,
            cex = 1,
            cex.axis = 0.7,
            cex.lab = 0.7
        )
        points(
            assessModRes$rank,
            assessModRes$modularity_max,
            xlab = "rank",
            ylab = "Q",
            pch = 1,
            cex = 1,
            cex.axis = 0.7,
            cex.lab = 0.7
        )
        abline(h=0.33, lty=2)
        
        par(new = TRUE)                             # Add new plot
        #par(mgp = (c(2, .1, 0)))
        plot(
            assessModRes$rank,
            assessModRes$n_community,
            pch = 17,
            cex = 1,
            axes = FALSE,
            xlab = "",
            ylab = ""
        )
        points(
            assessModRes$rank,
            assessModRes$n_community_max,
            pch = 2,
            cex = 1,
            xlab = "",
            ylab = ""
        )
        
        axis(side = 4, at = pretty(range(assessModRes$n_community)), cex.axis = 0.7)      # Add second axis
        mtext("n", side = 4, line = 1, cex=0.7)
        
        par(mar = (c(4, .1, 0.1, 0.1)))
        plot.new()
        legend(
            "bottom",
            c("Q", "Qmax", "n", "nmax"),
            col = "black",
            pch = c(16, 17, 1, 2),
            cex = 0.7
        )
        
    }