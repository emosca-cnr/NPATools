#' Plot modularity trend
#' @param assessModRes result of `assess_modularity()`
#' @param out_dir output directory
#' @return plot
#' @export
#' @importFrom graphics layout par points plot.new legend abline mtext
#' @importFrom grDevices png dev.off

plot_modu_trend <- function(assessModRes = NULL, out_dir=NULL) {
        
        if(is.null(out_dir)){
                out_dir <- getwd()
        }else{
                dir.create(out_dir, recursive = T, showWarnings = F)
        }
        
        png(file.path(out_dir, "mod_trend.png"), width = 180, height = 90, units="mm", res=300)
        
        #layout(matrix(1:3, ncol = 3), widths = c(0.45, .45, 0.1))
        par(mfrow = c(1, 2))
        par(mar = c(2.5, 2.5, 1, 2))
        par(mgp = (c(1.5, .5, 0)))
        
        plot(
            assessModRes$rank,
            assessModRes$modularity,
            xlab = "rank",
            ylab = "Q",
            pch = 16,
            cex = 0.7,
            cex.axis = 0.7,
            cex.lab = 0.7,
            lty=2,
            type="b"
        )
        points(
            assessModRes$rank,
            assessModRes$modularity_max,
            xlab = "rank",
            ylab = "Q",
            pch = 1,
            cex = 0.7,
            cex.axis = 0.7,
            cex.lab = 0.7,
            lty=2,
            type="b"
        )
        abline(h=0.33, lty=2)
        
        legend(
                "bottomright",
                c("all", "largest"),
                col = "black",
                pch = c(16, 1),
                cex = 0.7
        )
        
        #par(new = TRUE)                             # Add new plot
        #par(mgp = (c(2, .1, 0)))
        plot(
            assessModRes$rank,
            assessModRes$n_community,
            pch = 16,
            cex = 0.7,
            #axes = FALSE,
            xlab = "rank",
            ylab = "n",
            lty=2,
            type="b"
        )
        points(
            assessModRes$rank,
            assessModRes$n_community_max,
            pch = 1,
            cex = 0.7,
            xlab = "",
            ylab = "",
            lty=2,
            type="b"
        )
        
        legend(
                "bottomright",
                c("all", "largest"),
                col = "black",
                pch = c(16, 1),
                cex = 0.7
        )
        dev.off()
        
        #axis(side = 4, at = pretty(range(assessModRes$n_community)), cex.axis = 0.7)      # Add second axis
        #mtext("n", side = 4, line = 1, cex=0.7)
        
        par(mar = (c(4, .1, 0.1, 0.1)))
        plot.new()
        legend(
            "bottom",
            c("Q", "Qmax", "n", "nmax"),
            col = "black",
            pch = c(16, 1, 17, 2),
            cex = 0.7
        )
        dev.off()
        
    }