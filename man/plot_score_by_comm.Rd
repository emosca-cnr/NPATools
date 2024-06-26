% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_score_by_comm.R
\name{plot_score_by_comm}
\alias{plot_score_by_comm}
\title{Boxplots of a score by community}
\usage{
plot_score_by_comm(
  g = NULL,
  comm_attr = "comm_id",
  score_attr = NULL,
  min_comm_size = 5,
  top_genes_to_label = 3
)
}
\arguments{
\item{g}{igraph object}

\item{comm_attr}{attrubute that contains community information}

\item{score_attr}{attrubute that contains gene scores; high values are considered significant}

\item{min_comm_size}{minimum community size}

\item{top_genes_to_label}{number of genes to label}
}
\description{
Boxplots of a score by community
}
