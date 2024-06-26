% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exprs_to_cdf.R
\name{exprs_to_cdf}
\alias{exprs_to_cdf}
\title{Gene expression to cdf values}
\usage{
exprs_to_cdf(gene_expr = NULL, type = c("kcdf", "ecdf"))
}
\arguments{
\item{gene_expr}{normalized gene expression matrix}

\item{type}{"kcdf" will use Gaussian kernel while "ecdf" will use empirical estimate}
}
\description{
Gene expression to cdf values
}
