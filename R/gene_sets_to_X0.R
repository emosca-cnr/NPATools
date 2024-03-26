#' Quantify the occurence of gene sets among positive element of X0
#' @param X0 X0 matrix
#' @param gsl gene set list. The element must match the row names of X0
#' @param add_names TRUE/FALSE, whether to add a string with collapsed rown names of X0; if not only numbers will be included in the output
#' @export
#' 

gene_sets_to_X0 <- function(X0=NULL, gsl=NULL, add_names=F){
  
  if(add_names){
    ans <- lapply(gsl, function(x) X0[rownames(X0) %in% x, ])
    ans <- lapply(ans, function(x) c(colSums(x), unlist(apply(x, 2, function(y) paste0(sort(rownames(x)[y != 0]), collapse = ";")))))
    ans <- do.call(rbind, ans)
    colnames(ans) <- c(paste0(colnames(ans)[1:ncol(X0)], "_n"), paste0(colnames(ans)[1:ncol(X0)], "_id"))
  }else{
    ans <- lapply(gsl, function(x) colSums(X0[rownames(X0) %in% x, ]))
    ans <- do.call(rbind, ans)
    colnames(ans) <- paste0(colnames(ans), "_n")
  }
  
  ans <- as.data.frame(ans)
  
  return(ans)
}