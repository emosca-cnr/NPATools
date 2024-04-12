#' empirical False Discovery Rate
#' @param real_values vactor of real (observed) values
#' @param all_values vactor of all values (real + permuted) values
#' @param BPPARAM deprecated but kept for compatibility
#' @param correct.max logical, to decide if to correct the maximum values to 1 or not
#' @import BiocParallel
#' @export
eFDR <- function(real_values=NULL, all_values=NULL, BPPARAM = NULL, correct.max = TRUE){
  
  #FDR = (# x > permuted values / # x > real values) * (#real values / #permuted values)
  #approch used in GSEA
  
  # if (is.null(BPPARAM)) {
  #   BPPARAM <- SerialParam()
  # }
  # 
  # cat("BPPARAM\n")
  # print(BPPARAM)
  
  names(real_values) <- 1:length(real_values) #to keep initial order
  
  r_a_ratio <- length(real_values) / length(all_values)
  
  #fdr_values <- r_a_ratio * unlist(bplapply(real_values, function(x) sum(all_values >= x) / sum(real_values >= x), BPPARAM = BPPARAM))
  
  #By definition, the rank of a decreasing list is equal to the number of values higher than each element of the list
  real_values <- sort(real_values, decreasing = T)
  all_values <- sort(all_values, decreasing = T)
  
  real_values_rank <- rank(-real_values, ties.method = "max") 
  
  which_eq_real <- match(real_values, all_values) #place of each real value in all_value
  all_values_rank <- rank(-all_values, ties.method = "max")
  
  fdr_values <- r_a_ratio *  all_values_rank[which_eq_real] / real_values_rank
  
  #to obtain the initial order
  fdr_values <- fdr_values[match(as.character(1:length(real_values)), names(real_values_rank))]
  
  if(correct.max) {
    fdr_values[fdr_values>1] <- 1
  }

  return(fdr_values)
  
}