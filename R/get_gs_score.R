#' Get a gene set score from Xs matrix
#' @param Xs Xs matrix
#' @param gsl named list of gene sets, which must match the row names of X
#' @param FUN function to use to calculate the gene set score from the element of Xs; mean() by default
#' @export

get_gs_score <- function(Xs=NULL, gsl=NULL, FUN=NULL){
  
  if(is.null(FUN)){
    FUN <- mean
  }
  
  FUN <- match.fun(FUN)
  
  ans <- do.call(rbind, lapply(gsl, function(gs) apply(Xs[rownames(Xs) %in% gs, , drop=FALSE], 2, FUN=FUN)))
  
  return(ans)
  
}