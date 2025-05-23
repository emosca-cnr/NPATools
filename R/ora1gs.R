#' Hypergeometric test on 1 dataset
#' @param wb white balls
#' @param bb black balls
#' @param bd balls drawn
#' @importFrom stats phyper dhyper
#'

ora1gs <- function(wb=NULL, bb=NULL, bd=NULL){
  wbd <- bd[bd %in% wb]

  p_out <- phyper(length(wbd), length(wb), length(bb), length(bd), lower.tail = FALSE) + dhyper(length(wbd), length(wb), length(bb), length(bd))

  return(data.frame(wb=length(wb), bb=length(bb), bd=length(bd), wbd=length(wbd), p_val=p_out, stringsAsFactors = F))
}
