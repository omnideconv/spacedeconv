#' Cluster spacedeconv results
#'
#' Performs clustering on data obtained from a SpatialExperiment object using
#' specified methods. This function allows clustering based on deconvolution
#' results, expression data, pathway, or transcription factors (TF) analyses.
#'
#' @param spe A SpatialExperiment object containing the data to be clustered.
#' @param method A character vector specifying the clustering method to use.
#'        Options are "kmeans" and "hclust". Default is c("kmeans", "hclust").
#' @param data A character vector indicating the type of data to cluster.
#'        Options are "deconvolution", "expression", "pathway", or "tf".
#'        Default is c("deconvolution", "expression", "pathway", "tf").
#' @param dist_method A character vector specifying the distance measure to be
#'        used if method = "hclust". Default options are "correlation",
#'        "euclidean", "maximum", "manhattan", "canberra", "binary",
#'        "minkowski".
#' @param hclust_method A character vector indicating the agglomeration method
#'        to be used if method = "hclust". Default options are "complete",
#'        "ward.D", "ward.D2", "single", "average", "mcquitty", "median",
#'        "centroid".
#' @param nclusters An integer specifying the number of clusters to create.
#'        Default is 3.
#' @param spmethod A character string indicating the spatial method used for
#'        the clustering. Must be one of "dorothea", "progeny", "expression",
#'        or the name of the deconvolution method used. Default is NULL.
#' @param pca_dim An integer vector specifying PCA dimensions to be used for
#'        clustering of expression data. This is used with Seurat::FindNeighbors.
#'        Default is seq(1, 30).
#' @param clusres A numeric value for the clustering resolution to be used
#'        with Seurat::FindClusters. Default is 0.5.
#' @param ... Additional parameters.
#'
#' @return The modified SpatialExperiment object with clustering results.
#' @export
#'
cluster <- function(spe,
                    method = c("kmeans", "hclust"),
                    data = c("deconvolution", "expression", "pathway", "tf"),
                    dist_method = c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                    hclust_method = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"),
                    nclusters = 3,
                    spmethod = NULL,
                    pca_dim = seq(1, 30),
                    clusres = 0.5, ...) {
  cli::cli_rule(left = "spacedeconv")
  cli::cli_progress_step("testing parameter", msg_done = "parameter OK")

  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }

  cli::cli_progress_step("Extracting data", msg_done = "Extracted data for clustering")

  # convert to sparse matrices
  spe <- check_datatype(spe)

  data <- match.arg(data)

  # if (!data %in% c("expression", "deconvolution", "pathway", "tf")) {
  #   stop("`data` must be one of the following: expression, deconvolution, pathway, tf")
  # }

  if (data == "expression") {
    # convert the spe to a seurat object
    # create expression data matrix and rename the rows to normal gene names instead of the ENSEMBL
    expression_data <- spe@assays@data@listData[["counts"]]

    # create spatialcoordinates matrix
    spatial_coordinates <- as.matrix(SpatialExperiment::spatialCoords(spe))

    cli::cli_alert_info(paste("Clustering:", data))
    cli::cli_alert_info(paste("Cluster resolution:", clusres))



    cli::cli_progress_step("Extracting data", msg_done = "Extracted data for clustering")

    # create seurat object
    seurat_obj <- SeuratObject::CreateSeuratObject(
      counts = expression_data,
      spatial = spatial_coordinates,
      project = "ST"#,
      #assay = "Spatial"
    )

    # normalize spatial counts
    seurat_obj <- Seurat::SCTransform(seurat_obj,
      #assay = "Spatial",
      verbose = FALSE
    )

    # run PCA and do the clustering
    seurat_obj <- Seurat::RunPCA(seurat_obj,
      assay = "SCT",
      verbose = FALSE
    )
    seurat_obj <- Seurat::FindNeighbors(seurat_obj,
      reduction = "pca",
      dims = pca_dim
    )
    # clusres to define the resolution of the clustering
    for (i in clusres) {
      seurat_obj <- Seurat::FindClusters(seurat_obj,
        verbose = FALSE,
        resolution = i
      )
      cluster <- seurat_obj@meta.data[["seurat_clusters"]]
      names(cluster) <- rownames(seurat_obj@meta.data)
      cname <- paste0("cluster_expression_", i)
      SummarizedExperiment::colData(spe)[cname] <- cluster
      cli::cli_progress_update()
    }
  } else if (data %in% c("deconvolution", "tf", "pathway")) {
    if (is.null(spmethod)) {
      if (data == "tf") {
        spmethod <- "dorothea"
      } else if (data == "pathway") {
        spmethod <- "progeny"
      } else {
        stop("Parameter 'spmethod' is null or missing, but is required")
      }
    }

    dist_method <- match.arg(dist_method)
    hclust_method <- match.arg(hclust_method)
    method <- match.arg(method)

    cli::cli_alert_info(paste("Clustering:", data))
    cli::cli_alert_info(paste("By:", method))
    cli::cli_alert_info(paste("Number of clusters:", nclusters))
    if (method == "hclust") {
      cli::cli_alert_info(paste("Distance Method:", dist_method))
      cli::cli_alert_info(paste("Hclust Method:", hclust_method))
    }

    tmp <- SummarizedExperiment::colData(spe)[available_results(spe,
      method = spmethod
    )]
    # clusters
    cli::cli_progress_bar("Clustering", total = length(nclusters))
    for (i in nclusters) {
      result <- switch(method,
        kmeans = {
          stats::kmeans(tmp, centers = i, ...)$cluster
        },
        hclust = {
          if (dist_method == "correlation") {
            d <- as.dist(1 - cor(t(as.matrix(tmp))))
          } else {
            d <- dist(tmp, method = dist_method)
          }
          hc <- stats::hclust(d, method = hclust_method)
          stats::cutree(hc, k = i)
        }
      )

      cluster <- as.factor(result)
      cname <- paste("cluster", spmethod, i, sep = "_")
      SummarizedExperiment::colData(spe)[cname] <- cluster
      cli::cli_progress_update()
    }
  }



  cli::cli_progress_done()

  return(spe)
}

#' Get Top Features
#' @param idx ids
#' @param scores scores
#' @param topn number of features
topfeat <- function(idx, scores, topn) {
  if (length(idx) > 1) {
    return(head(sort(apply(scores[idx, ], 2, mean), decreasing = TRUE), topn))
  } else {
    return(head(sort(scores[idx, ], decreasing = TRUE), topn))
  }
}

#' Get cluster features
#' @param spe spatialExperiment with cluster results
#' @param clusterid = name of the column with the clustering results
#' @param topn number of top features to be shown
#' @param spmethod spatial method used fot the clustering, must be dorothea, progeny, expression or the name of the deconvolution method used
#' @param zscore = if the results should be z-score scaled or not
#' @export
get_cluster_features <- function(spe,
                                 clusterid = NULL,
                                 topn = 3,
                                 spmethod = NULL,
                                 zscore = TRUE) {
  # Many checks are still missing
  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }
  if (is.null(clusterid)) {
    stop("Parameter 'clusterid' is null or missing, but is required")
  }
  if (is.null(spmethod)) {
    spmethod <- unlist(strsplit(clusterid, "_"))[2]
  }

  # check if spmethod is actually available
  if (!any(grepl(spmethod, names(colData(spe))))) {
    stop(paste("spmethod", spmethod, "not found in spe object"))
  }

  # Extract clusters
  clusters <- colData(spe)[, available_results(spe, method = "cluster"),
    drop = FALSE
  ]
  clusters <- clusters[, clusterid]



  # Scores
  if (spmethod == "expression") {
    if (is.element("cpm", assayNames(spe)) && !is.null(spe@assays@data$cpm)) {
      # if all good, extract the CPM values from the spe in a new expression_data variable
      scores <- t(as.matrix(spe@assays@data$cpm))
    } else {
      # normalize the spe object and then extract the cpm
      spe <- spacedeconv::normalize(spe, method = "cpm", assay = "counts")

      # extract the cpm as above in the scores variable (take as.matrix to avoid sparse matrix)
      scores <- t(as.matrix(spe@assays@data$cpm))
    }
  } else {
    scores <- colData(spe)[, available_results(spe, method = spmethod)]
  }

  # Transform to z-scores
  if (zscore) {
    scores <- as.matrix(scores)

    # number of complete cases
    num_na_rows <- sum(!complete.cases(scores))

    # remove incomplete cases from scores and clusters
    clusters <- clusters[complete.cases(scores)]
    scores <- scores[complete.cases(scores), ]

    # notify user
    if (num_na_rows > 0) {
      cli::cli_alert_info(paste(num_na_rows, "rows were removed due to NA values."))
    }

    scores <- t((t(scores) - apply(scores, 2, mean)) / apply(scores, 2, sd))
  }

  # Compute most abundant cell types
  idx <- split(names(clusters), clusters)
  topscores <- lapply(idx, topfeat, scores = scores, topn = topn)

  return(topscores)
}
