% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decoupleR.R
\name{compute_activities}
\alias{compute_activities}
\title{Compute Pathway and Transcription Factor Activities and integrate into SpatialExperiment Using decoupleR}
\usage{
compute_activities(
  spe,
  reference,
  method = "wmean",
  assay = "cpm",
  statistic = NULL,
  ...
)
}
\arguments{
\item{spe}{The \code{SpatialExperiment} object containing spatially resolved expression data.}

\item{reference}{A reference dataframe, computed using \code{get_decoupleR_reference()}}

\item{method}{The method used to compute activities, such as "wmean", "viper", "gsva", etc. see \code{decoupleR} documentation for details.}

\item{assay}{The name of the assay within the \code{SpatialExperiment} to use for calculations.}

\item{statistic}{In cases where the chosen method produces multiple results, this parameter allows selection of a specific statistic for inclusion in the output.}

\item{...}{Further arguments passed to the computational methods.}
}
\value{
The function returns the input \code{SpatialExperiment} object with additional metadata columns representing the computed activities.
}
\description{
This function computes pathway activities or transcription factor activities using a \code{SpatialExperiment} object
using various computational methods provided by the \code{decoupleR} package.
}
