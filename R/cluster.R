#' Cluster Spots
#'
#' Cluster spots and add cluster annotation to the spatial object
#'
#' @param spe SpatialExperiment
#' @param method clustering method, (kmeans)
#' @param assay assay to use for clustering
#' @param nclusters number of clusters, for kmeans
#' @param ... further arguments passed to the clustering functions
#'
#' @returns SpatialExperiment containing cluster information
#'
cluster <- function(spe, method="kmeans", assay = "counts", nclusters = 3, ...){
  if (is.null(spe)){
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if (!assay %in% names(SummarizedExperiment::assays(spe))){
    stop("Requested assay not available in provided SpatialExperiment")
  }

  # assay
  tmp <- SummarizedExperiment::assay(spe, assay)
  tmp <- t(as.matrix(tmp)) # cluster spots and not genes

  ########### remove as.matrix!!!!!

  # cluster
  cluster <- switch (method,
    kmeans = {stats::kmeans(tmp, centers = nclusters, ...)$cluster}
  )

  # add clustering
  cluster <- as.factor(cluster)
  colData(spe) <- cbind(colData(spe), cluster)

  return (spe)
}
