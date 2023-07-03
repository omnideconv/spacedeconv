#' Cluster spacedeconv results
#' @param spe SpatialExperiment
#' @param method clustering method to be chosen between kmeans and hclust
#' @param data what data to cluster
#' @param dist_method the distance measure to be used if method = hclust
#' @param hclust_method the agglomeration method to be used if method = hclust
#' @param nclusters number of clusters
#' @param spmethod spatial method used fot the clustering, must be dorothea, progeny, expression or the name of the deconvolution method used
#' @param pca_dim PCA dimensions to be used for the clustering of expression data - Seurat::FindNeighbors
#' @param clusres clustering resolution to be used for the clustering of expression data - Seurat::FindClusters


#' @export

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

  if (data == "expression") {
    # convert the spe to a seurat object
    # create expression data matrix and rename the rows to normal gene names instead of the ENSEMBL
    expression_data <- spe@assays@data@listData[["counts"]]

    # create spatialcoordinates matrix
    spatial_coordinates <- as.matrix(SpatialExperiment::spatialCoords(spe))



    cli::cli_progress_step("Extracting data", msg_done = "Extracted data for clustering")

    # create seurat object
    seurat_obj <- SeuratObject::CreateSeuratObject(
      counts = expression_data,
      spatial = spatial_coordinates,
      project = "ST",
      assay = "Spatial"
    )

    # normalize spatial counts
    seurat_obj <- Seurat::SCTransform(seurat_obj,
      assay = "Spatial",
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
      stop("Parameter 'spmethod' is null or missing, but is required")
    }

    if (!data %in% c("expression", "deconvolution", "pathway", "tf")) {
      stop("`data` must be one of the following: expression, deconvolution, pathway, tf")
    }

    if (!dist_method %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
      stop("`dist_method` must be one of the following: correlation, euclidean, maximum, manhattan, canberra, binary or minkowski")
    }

    if (!hclust_method %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")) {
      stop("`hclust_method` must be one of the following: ward.D, ward.D2, single, complete, average, mcquitty, median,centroid")
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
    scores <- t((t(scores) - apply(scores, 2, mean)) / apply(scores, 2, sd))
  }

  # Compute most abundant cell types
  idx <- split(names(clusters), clusters)
  topscores <- lapply(idx, topfeat, scores = scores, topn = topn)

  return(topscores)
}
