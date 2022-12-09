#' Cluster Spots
#'
#' Cluster spots and add cluster annotation to the spatial object
#'
#' @param spe SpatialExperiment
#' @param method clustering method, (kmeans)
#' @param cluster what to cluster, c("expression", "deconvolution")
#' @param assay assay to use for clustering
#' @param nclusters number of clusters, for kmeans
#' @param ... further arguments passed to the clustering functions
#'
#' @returns SpatialExperiment containing cluster information
#'
cluster <- function(spe, method="kmeans", cluster = "expression", assay = "counts", nclusters = 3, ...){
  if (is.null(spe)){
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if (!assay %in% names(SummarizedExperiment::assays(spe))){
    stop("Requested assay not available in provided SpatialExperiment")
  }

  if (!cluster %in% c("expression", "deconvolution")){
    stop("can only cluster expression or deconvolution")
  }

  if (cluster == "expression"){
    # assay
    tmp <- SummarizedExperiment::assay(spe, assay)
    tmp <- t(as.matrix(tmp)) # cluster spots and not genes

    ##### remove as.matrix!!!!!

  } else if (cluster == "deconvolution"){
    tmp <- SummarizedExperiment::colData(spe)[available_results(spe)]
  }

  # cluster
  cluster <- switch (method,
    kmeans = {stats::kmeans(tmp, centers = nclusters, ...)$cluster}
  )

  # add clustering
  cluster <- as.factor(cluster)
  SummarizedExperiment::colData(spe) <- cbind(SummarizedExperiment::colData(spe), cluster)

  return (spe)
}
