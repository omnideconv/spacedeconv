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
#' @export
#'
cluster <- function(spe, method = "kmeans", cluster = "expression", assay = "counts", nclusters = 3, ...) {
  cli::cli_rule(left="spacedeconv")

  cli::cli_progress_step("testing parameter", msg_done = "parameter OK")

  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if (!assay %in% names(SummarizedExperiment::assays(spe))) {
    stop("Requested assay not available in provided SpatialExperiment")
  }

  if (!cluster %in% c("expression", "deconvolution")) {
    stop("can only cluster expression or deconvolution")
  }

  cli::cli_progress_step("Extracting data", msg_done = "Extracted data for clustering")

  if (cluster == "expression") {
    # assay
    tmp <- SummarizedExperiment::assay(spe, assay)
    tmp <- t(as.matrix(tmp)) # cluster spots and not genes

    ##### remove as.matrix!!!!!
  } else if (cluster == "deconvolution") {
    tmp <- SummarizedExperiment::colData(spe)[available_results(spe)]
  }

  cli::cli_progress_step("Clustering", msg_done = "Finished clustering")

  # cluster
  cluster <- switch(method,
    kmeans = {
      stats::kmeans(tmp, centers = nclusters, ...)$cluster
    }
  )

  # add clustering
  cluster <- as.factor(cluster)
  SummarizedExperiment::colData(spe) <- cbind(SummarizedExperiment::colData(spe), cluster)

  cli::cli_progress_done()

  return(spe)
}

#' Get Most Abundant cell types per cluster
#'
#' @param spe SpatialExperiment with deconvolution results
#' @param method how to calculate (sum, mean)
#' @param k return top k most abundant cell types
#' @export
get_mostAbundantInCluster <- function(spe, method="mean", k=3){
  df <- colData(spe)[, available_results(spe)]
  clusters <- df$cluster
  # df$cluster <- NULL

  if (!method %in% c("mean", "median")){
    stop("Method not supported")
  }

  result <- list()

  for (cluster in unique(clusters)){
    tmp <- df[df$cluster==cluster, ]
    tmp$cluster <- NULL

    res <- sort(apply(tmp, 2, method), decreasing = T)
    res <- res[1:k]
    result[[cluster]] <- res
  }

  return (result)
}
