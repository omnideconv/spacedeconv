% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster.R
\name{cluster}
\alias{cluster}
\title{Cluster spacedeconv results}
\usage{
cluster(
  spe,
  method = c("kmeans", "hclust"),
  spmethod = c("expression", "progeny", "dorothea", "collectri",
    unname(deconvolution_methods)),
  dist_method = c("correlation", "euclidean", "maximum", "manhattan", "canberra",
    "binary", "minkowski"),
  hclust_method = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty",
    "median", "centroid"),
  nclusters = 3,
  pca_dim = seq(1, 30),
  clusres = 0.5,
  ...
)
}
\arguments{
\item{spe}{A SpatialExperiment object containing the data to be clustered.}

\item{method}{A character vector specifying the clustering method to use.
Options are "kmeans" and "hclust". Default is c("kmeans", "hclust").}

\item{spmethod}{A character vector indicating the type of analysis for
clustering. Options include "expression", "progeny", "dorothea",
"collectri", and names of deconvolution methods. Default is based on
the available deconvolution methods.}

\item{dist_method}{A character vector specifying the distance measure to be
used for "hclust" method. Options include "correlation", "euclidean",
"maximum", "manhattan", "canberra", "binary", "minkowski".}

\item{hclust_method}{A character vector indicating the agglomeration method
to be used with "hclust". Options include "complete", "ward.D",
"ward.D2", "single", "average", "mcquitty", "median", "centroid".}

\item{nclusters}{An integer specifying the number of clusters to create.
Default is 3.}

\item{pca_dim}{An integer vector specifying PCA dimensions to be used for
clustering of expression data. This is used with Seurat::FindNeighbors.
Default is seq(1, 30).}

\item{clusres}{A numeric value for the clustering resolution to be used
with Seurat::FindClusters. Default is 0.5.}

\item{...}{Additional parameters for clustering methods.}
}
\value{
A modified SpatialExperiment object with added clustering results.
}
\description{
Performs clustering on data obtained from a SpatialExperiment object using
specified methods. This function allows for clustering based on deconvolution
results, expression data, pathway, or transcription factors (TF) analyses.
}
