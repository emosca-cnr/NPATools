#' Select the most representative features
#' @param gsl gene set list
#' @param method on of the method available for proxy::simil()
#' @importFrom proxy simil
#' @export

calc_gs_sim <- function(gsl=NULL, method="Simpson"){
  
  feat_occ <- feat_sets_to_occ(gsl)
  oi <- simil(as.matrix(feat_occ), method = method, diag = F, upper = F)
  
  return(oi)

}
