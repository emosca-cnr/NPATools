% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comm_net.R
\name{comm_net}
\alias{comm_net}
\title{Create a network of communities}
\usage{
comm_net(g = NULL, remove.multiple = FALSE)
}
\arguments{
\item{g}{igraph object, vertex attribute "comm_id" is mandatory}

\item{remove.multiple}{whether to remove multiple links among communities}
}
\value{
igraph object
}
\description{
Function for reduce the network as a community network
}
