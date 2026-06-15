#' Provide a summary score for he networks composed of a ranked list of genes
#' @param net_scores the table provided by the function score_networks()
#' @param LCC whether to consider the largest connected components (TRUE) instead of considering all connected components (FALSE, default)
#' @param minpeakdistance the minimum distance (in indices) peaks have to have to be counted
#' @param minpeakheight the minimum (absolute) height a peak has to have to be recognized as such
#' @param ... further argumnts for pracma::findpeaks()
#' @returns A data frame with:
#' * the input scores normalized to their maximum; `not_linkers` is 1 minus linkers (divided to their maximum);
#' * `y`: the summary score, which is the sum of the previous columns;
#' * `peak`: the presence of a peak in y.
#' @export
#' @importFrom stats setNames
#' @importFrom pracma findpeaks

score_networks_summary <- function(net_scores = NULL, LCC=FALSE, minpeakdistance=1, minpeakheight=0, ...) {
	
	stopifnot(!is.null(net_scores))
	
	
	ans <- data.frame(rank=net_scores$rank)
	
	ans$z_n <- net_scores$z / max(net_scores$z)
	
	if(LCC){
		
		ans$X0_n <- net_scores$X0_LCC / max(net_scores$X0_LCC)
		ans$Q_n <- net_scores$Q_LCC / max(net_scores$Q_LCC)
		
	}else{
		
		ans$X0_n <- net_scores$X0_CC / max(net_scores$X0_CC)
		ans$Q_n <- net_scores$Q_CC / max(net_scores$Q_CC)
		
	}
	
	ans$not_linkers <- 1 - net_scores$linkers / max(net_scores$linkers)

	ans$y <- ans$z_n + ans$X0_n + ans$Q_n + ans$not_linkers
	
	peaks <- as.data.frame(findpeaks(ans$y, minpeakdistance = minpeakdistance, minpeakheight=minpeakheight, ...))
	colnames(peaks) <- c("height", "idx", "from", "to")
	peaks$peak <- TRUE
	peaks$rank <- ans$rank[peaks$idx]
	peaks$height <- peaks$idx <- peaks$from <- peaks$to <- NULL
	
	ans <- merge(ans, peaks, by="rank", all.x=T)
	ans$peak[is.na(ans$peak)] <- FALSE
	if(!ans$peak[which.max(ans$y)]){ #if the last point is the highest findpeaks does not consider it as a peak
		ans$peak[which.max(ans$y)] <- TRUE
	}
	
	ans <- ans[order(as.numeric(ans$rank)), ]
	
	return(ans)
	
}

