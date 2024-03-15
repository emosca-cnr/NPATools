#' Normalization of Nd values
#' @export
normalize_nd <- function(x=NULL, na.rm = TRUE, method=c("max", "min_max", "sum"), ymax=1, ymin=0) {
  
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
