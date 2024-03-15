#' Gene expression to cdf values
#' @export
#' @importFrom spatstat.explore CDF
#' @importFrom stats ecdf

exprs_to_cdf <- function(gene_expr=NULL, type=c("kcdf", "ecdf")){
  
  col_names <- colnames(gene_expr)
  gene_expr <- t(apply(gene_expr, 1, get_cdf, type=type))
  colnames(gene_expr) <- col_names
  print(dim(gene_expr))
  
  print(summary(apply(gene_expr, 1, min)))
  print(summary(apply(gene_expr, 1, max)))

  
  return(gene_expr)
  
}