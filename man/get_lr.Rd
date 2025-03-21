% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_lr.R
\name{get_lr}
\alias{get_lr}
\title{Compute L-R score for each spot}
\usage{
get_lr(spe, resource = "Consensus", method = "min", organism = "human")
}
\arguments{
\item{spe}{SpatialExperiment object}

\item{resource}{"consensus" table from Omnipath as default, the user can provide a data frame containing L-R pairs. The data frame should contain at least the following two columns: source_genesymbol = ligands, target_genesymbol = receptors}

\item{method}{mathematical approach to compute L-R scores. Options: min, and product L-R, default is product}

\item{organism}{choose the organism to be considered, default human, options: human or mouse}
}
\description{
Compute L-R score for each spot
}
