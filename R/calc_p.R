#' Estimation of p values
#' @description This function is used to calculate p-value from a list of matrices
#' @details the functions calculate p-value for a permutation-based approach. It takes as an input a list of matrices, where the first is the one with the real data
#' and the other are the data obtained through permutations. Then, the function searches for how many values in the permutations are equal or higher to the real.
#' These values are then divided by the number matrices (1 + number of permutations)
#' @param X list of matrices, where the first is the one obtained with real data and the others are obtained through permutations. 
#' The matrices should be 1 column matrices with the same order.
#' @param type g for >=, l <= than the real values
#' @export
#'
calc_p <- function(X=NULL, type=c("g", "l")){
  
  type <- match.arg(type, c("g", "l"))
  stopifnot(length(X)>1)
  
  p <- matrix(0, nrow=nrow(X[[1]]), ncol=ncol(X[[1]]), 
              dimnames = list(rownames(X[[1]]), colnames(X[[1]])))
  
  if(type=="g"){
    p <- lapply(X, function(x) x >= X[[1]])
  }
  if(type=="l"){
    p <- lapply(X, function(x) x <= X[[1]])
  }
  p <- Reduce("+", p)

  p <- p / length(X)
  
  return(p)
}