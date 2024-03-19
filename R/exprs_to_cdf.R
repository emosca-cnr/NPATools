#' Gene expression to cdf values
#' @param gene_expr normalized gene expression matrix
#' @param type "kcdf" will use Gaussian kernel while "ecdf" will use empirical estimate
#' @export

exprs_to_cdf <- function(gene_expr=NULL, type=c("kcdf", "ecdf")){
  
  col_names <- colnames(gene_expr)
  gene_expr <- t(apply(gene_expr, 1, get_cdf, type=type))
  colnames(gene_expr) <- col_names
  print(dim(gene_expr))
  
  print(summary(apply(gene_expr, 1, min)))
  print(summary(apply(gene_expr, 1, max)))

  return(gene_expr)
  
}