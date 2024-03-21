#' empirical False Discovery Rate
#' @param real_values vactor of real (observed) values
#' @param all_values vactor of all values (real + permuted) values
#' @param BPPARAM An optional BiocParallelParam instance determining the parallel back-end to be used during evaluation. If NULL, parallel evaluation is disabled using SerialParam(). See ?bplapply.
#' @param correct.max logical, to decide if to correct the maximum values to 1 or not
#' @import BiocParallel
#' @export
eFDR <- function(real_values=NULL, all_values=NULL, BPPARAM = NULL, correct.max = TRUE){
  
  #FDR = (# x > permuted values / # x > real values) * (#real values / #permuted values)
  #approch used in GSEA
  
  if (is.null(BPPARAM)) {
    BPPARAM <- SerialParam()
  }
  
  cat("BPPARAM\n")
  print(BPPARAM)
  
  fdr_values <- length(real_values) / length(all_values)
  fdr_values <- fdr_values * unlist(bplapply(real_values, function(x) sum(all_values >= x) / sum(real_values >= x), BPPARAM = BPPARAM))
  if(correct.max) {
    fdr_values[fdr_values>1] <- 1
  }

  return(fdr_values)
  
}