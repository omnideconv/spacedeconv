% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization.R
\name{normalize}
\alias{normalize}
\title{Normalize Gene Expression Data in SingleCellExperiment or SpatialExperiment Objects}
\usage{
normalize(object, method = "cpm", assay = "counts")
}
\arguments{
\item{object}{A \code{SingleCellExperiment} or \code{SpatialExperiment} object containing raw gene expression data.}

\item{method}{The normalization method to apply: "cpm" for counts per million, "logcpm" for log-transformed counts per million, or "logpf" for log-proportional factor.}

\item{assay}{The name of the assay to normalize, default "counts".}
}
\description{
This function applies normalization to the gene expression data stored within a \code{SingleCellExperiment} or \code{SpatialExperiment} object.
It supports both counts per million (CPM) and log-transformed counts per million (logCPM) normalization methods.
}
