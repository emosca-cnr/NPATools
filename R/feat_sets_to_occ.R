#' Feature set list to occurrence matrix
#' @param x gene set list
#' @importFrom Matrix Matrix
#' @export

feat_sets_to_occ <- function(x=NULL){
  
  all_x <- sort(unique(unlist(x)))
  
  ans <- Matrix(0, nrow = length(x), ncol = length(all_x), dimnames = list(names(x), all_x))
  
  for(i in 1:length(x)){
    ans[i, colnames(ans) %in% x[[i]]] <- 1
  }
  
  return(ans)
  
  
}