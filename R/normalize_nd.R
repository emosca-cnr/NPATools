#' Normalization of ND steady state values
#' @param x numeric vector with steady state diffusion values
#' @param method normalization approach
#' @param ymin min value for min-max normalization
#' @param ymax max value for min-max normalization
#' @return Normalized ND values
#' @export
normalize_nd <- function(x=NULL, method=c("max", "min_max", "sum"), ymax=1, ymin=0) {
  
  method <- match.arg(method, c("max", "min_max", "sum"))
  ans <- x
  
  if(any(x!=0)){
    
    if(method=="min_max"){
      ans <- (x- min(x)) /(max(x)-min(x)) * (ymax - ymin) + ymin
    }
    
    if(method=="max"){
      ans <- x/max(x)
    }
    
    if(method=="sum"){
      ans <- x/sum(x)
    }
  }
  return(ans)
}
