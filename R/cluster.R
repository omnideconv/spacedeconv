#' Cluster SpatialExperiment Results
#'
#' Performs clustering on a `SpatialExperiment` using either expression data or
#' selected result columns (deconvolution, pathway, or TF activity). Results are
#' stored as new `colData` columns.
#'
#' @param spe A `SpatialExperiment` object containing the data to be clustered.
#' @param method Clustering method: `"kmeans"` or `"hclust"`.
#' @param spmethod Data to cluster: `"expression"`, `"progeny"`, `"dorothea"`,
#' `"collectri"`, or a deconvolution method token.
#' @param dist_method Distance metric for `"hclust"`: `"correlation"`,
#' `"euclidean"`, `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"`,
#' `"minkowski"`.
#' @param hclust_method Agglomeration method for `"hclust"`: `"complete"`,
#' `"ward.D"`, `"ward.D2"`, `"single"`, `"average"`, `"mcquitty"`, `"median"`,
#' `"centroid"`.
#' @param nclusters Number of clusters to create.
#' @param pca_dim PCA dimensions used for expression clustering.
#' @param clusres Resolution for Seurat-based expression clustering.
#' @param ... Additional parameters for clustering methods.
#'
#' @return `SpatialExperiment` with added clustering results in `colData`.
#' @export
cluster <- function(spe,
                    method = c("kmeans", "hclust"),
                    spmethod = c("expression", "progeny", "dorothea", "collectri", unname(deconvolution_methods)),
                    dist_method = c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                    hclust_method = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid"),
                    nclusters = 3,
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

  spmethod <- match.arg(spmethod)

  # if (!data %in% c("expression", "deconvolution", "pathway", "tf")) {
  #   stop("`data` must be one of the following: expression, deconvolution, pathway, tf")
  # }

  if (spmethod == "expression") {
    # convert the spe to a seurat object
    # create expression data matrix and rename the rows to normal gene names instead of the ENSEMBL
    expression_data <- spe@assays@data@listData[["counts"]]

    # create spatialcoordinates matrix
    spatial_coordinates <- as.matrix(SpatialExperiment::spatialCoords(spe))

    cli::cli_alert_info(paste("Clustering:", spmethod))
    cli::cli_alert_info(paste("Cluster resolution:", toString(clusres)))


    cli::cli_progress_step("Extracting data", msg_done = "Extracted data for clustering")

    # create seurat object
    seurat_obj <- SeuratObject::CreateSeuratObject(
      counts = expression_data,
      ########## spatial = spatial_coordinates, # not using!
      project = "ST" # ,
      # assay = "Spatial"
    )

    # normalize spatial counts
    seurat_obj <- Seurat::SCTransform(seurat_obj,
      # assay = "Spatial",
      verbose = F
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
      cname <- paste0("cluster_expression_res_", i)
      SummarizedExperiment::colData(spe)[cname] <- cluster
      cli::cli_progress_update()
    }
  } else if (spmethod %in% c("progeny", "dorothea", "collectri", unname(deconvolution_methods))) {
    dist_method <- match.arg(dist_method)
    hclust_method <- match.arg(hclust_method)
    method <- match.arg(method)

    cli::cli_alert_info(paste("Clustering:", spmethod))
    cli::cli_alert_info(paste("By:", method))
    cli::cli_alert_info(paste("Number of clusters:", toString(nclusters)))
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
      cname <- paste("cluster", spmethod, "nclusters", i, sep = "_")
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

#' Get Cluster Features
#'
#' Summarizes top features per cluster from expression data or selected result
#' columns, with optional z-score scaling. Features are ranked by their mean
#' (or z-scored mean) within each cluster, and the top `topn` are returned.
#'
#' @param spe `SpatialExperiment` with clustering results.
#' @param clusterid Name of the column containing cluster labels.
#' @param topn Number of top features to return per cluster.
#' @param spmethod Method used for clustering (e.g., `dorothea`, `collectri`,
#' `progeny`, `expression`, or a deconvolution method token).
#' @param zscore Logical; z-score scale features before ranking.
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
  # clusters <- colData(spe)[, available_results(spe, method = "cluster"),
  #   drop = FALSE
  # ]
  # clusters <- clusters[, clusterid]
  clusters <- colData(spe)[, clusterid]


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
