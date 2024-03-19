#' Gene-level cdf estimation
#' @param gene_expr normalized gene expression profile
#' @param type "kcdf" will use Gaussian kernel while "ecdf" will use empirical estimate
#' @export
#' @importFrom spatstat.explore CDF
#' @importFrom stats ecdf density

get_cdf <- function(gene_expr=NULL, type=c("ecdf", "kcdf")){
  
  type <- match.arg(type, c("ecdf", "kcdf"))
  
  gene_expr <- as.numeric(gene_expr)
  
  if(type=="ecdf"){
    #empirical
    f_ecdf <- ecdf(gene_expr)
    ans <- f_ecdf(gene_expr)
    
  }else{
    #gaussian
    d <- density(gene_expr)
    f_kcdf <- CDF(d)
    ans <- f_kcdf(gene_expr)
  }
  
  return(ans)  

}
