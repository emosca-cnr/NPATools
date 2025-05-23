% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_net_enrich.R
\name{plot_net_enrich}
\alias{plot_net_enrich}
\title{Plot network enrichment results}
\usage{
plot_net_enrich(netEnrRes = NULL, showTopSig = 10, sigStat = "p_val")
}
\arguments{
\item{netEnrRes}{output of `assess_enrichment()`}

\item{showTopSig}{if not NULL, an integer that indicatesd the number of top significant ranks to display}

\item{sigStat}{the quantity to be plotted: for ORA, can be any of "er", "p" or "p_adj"; for GSEA can be "es", "p", "p_val", "nes" or "FDRq"}
}
\description{
Plot network enrichment results
}
