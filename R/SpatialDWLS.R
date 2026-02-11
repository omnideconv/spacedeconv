#' Build a SpatialDWLS Signature
#'
#' Constructs a DWLS signature matrix from single-cell data using Giotto and a
#' marker-selection strategy: either compute markers with a supported method
#' (`scran`, `gini`, `mast`) or pass a custom vector of marker genes directly.
#'
#' @param single_cell_obj `SingleCellExperiment` with cell type annotations.
#' @param assay_sc Single-cell assay to use.
#' @param marker_method Marker genes or method to compute markers (`scran`,
#' `gini`, or `mast`).
#' @param topNgenes Number of top-ranked marker genes per cluster.
#' @param cell_type_col Column containing cell type labels.
#' @param dim_method Dimension reduction method for Giotto workflow:
#' `pca`, `tsne`, or `umap`.
#' @param cluster_method Clustering method for Giotto workflow:
#' `leiden`, `louvain`, `kmeans`, or `hclust`.
#' @param ... Additional parameters passed to Giotto helpers.
build_model_spatial_dwls <- function(single_cell_obj, assay_sc = "counts", marker_method = "scran", topNgenes = 100, cell_type_col = "cell_ontology_class", dim_method = "pca", cluster_method = "leiden", ...) {
  # TODO Checks
  if (!checkCol(single_cell_obj, cell_type_col)) {
    stop("cell_type_col not available")
  }

  if (!exists("giotto_instructions")) {
    envname <- getOption("spacedeconv.conda_env", default = "spacedeconv-env")
    conda_envs <- reticulate::conda_list()
    python_path <- conda_envs$python[match(envname, conda_envs$name)]
    if (is.na(python_path)) {
      stop("Conda environment '", envname, "' not found.", call. = FALSE)
    }
    giotto_instructions <<- Giotto::createGiottoInstructions(python_path = python_path)
  }

  # check if requested assay exists
  if (!assay_sc %in% names(SummarizedExperiment::assays(single_cell_obj))) {
    message(
      "requested assay ",
      assay_sc,
      " not available in expression object. Using first available assay."
    )
    assay_sc <- names(SummarizedExperiment::assays(single_cell_obj))[1] # change to first available assay request not available
  }

  # spatial expression
  scExpression <- SummarizedExperiment::assay(single_cell_obj, assay_sc) %>% as("dgCMatrix")


  # markers
  if (marker_method %in% c("scran", "gini", "mast")) {
    message("Calculating markers by using the method: ", marker_method)
    obj <- Giotto::createGiottoObject(scExpression, instructions = giotto_instructions, ...)

    obj <- doGiottoWorkflow(obj, dim_method = dim_method, cluster_method = cluster_method)


    message("Calculating Markers")

    markers <- Giotto::findMarkers(obj, method = marker_method, cluster_column = paste0(cluster_method, "_clus"))

    # extract the top genes for each cluster
    genes <- character()
    for (cluster in markers) {
      genes <- c(genes, cluster$feats[1:topNgenes])
    }

    genes <- unique(genes)

    message("Giotto: extracted a total of ", length(genes), " marker genes")
  } else {
    message("Using the provided marker genes")
    genes <- marker_method
  }

  # get cell type vector from object
  cell_types <- as.vector(single_cell_obj[[cell_type_col]])

  signature <- Giotto::makeSignMatrixDWLSfromMatrix(matrix = scExpression, sign_gene = genes, cell_type_vector = cell_types, ...)

  return(signature)
}

#' Deconvolute Spatial DWLS
#'
#' @param spatial_obj Spatial Experiment
#' @param signature Signature
#' @param assay_sp Assay of SpatialExperiment to use
#' @param result_name token to identify deconvolution results in object, default = "spatialdwls"
#' @param ... additional parameters
deconvolute_spatial_dwls <- function(spatial_obj, signature, assay_sp = "counts", result_name = "spatialdwls", ...) {
  # create Giotto Object
  spExpression <- SummarizedExperiment::assay(spatial_obj, assay_sp) %>% as("dgCMatrix")
  spCoords <- SpatialExperiment::spatialCoords(spatial_obj)

  if (!exists("giotto_instructions")) {
    envname <- getOption("spacedeconv.conda_env", default = "spacedeconv-env")
    conda_envs <- reticulate::conda_list()
    python_path <- conda_envs$python[match(envname, conda_envs$name)]
    if (is.na(python_path)) {
      stop("Conda environment '", envname, "' not found.", call. = FALSE)
    }
    giotto_instructions <<- Giotto::createGiottoInstructions(python_path = python_path)
  }


  obj <- Giotto::createGiottoObject(expression = spExpression, spatial_locs = spCoords, instructions = giotto_instructions)

  obj <- doGiottoWorkflow(obj, ...)

  deconvolution <- Giotto::runDWLSDeconv(obj, sign_matrix = signature, n_cell = 10, return_gobject = TRUE)

  df <- data.frame(Giotto::getSpatialEnrichment(deconvolution, output = "data.table"))

  df$cell_ID <- NULL

  # attach method token
  df <- attachToken(df, result_name)

  return(df)
}

#' Perform Giotto Workflow of Normalization, HVG, Dimension Reduction and Clustering
#'
#' @param obj Giotto Object
#' @param calculateHVG Wheter to calculate HVG
#' @param dim_method dimension reduction method to use: pca, tsne, umap
#' @param cluster_method cluster method to use c(leiden, kmeans, hclust, louvain)
#' @param show_plot show plot
#' @param return_plot return plot
#' @param save_plot save plot
#' @param ... additional parameters
doGiottoWorkflow <- function(obj, calculateHVG = TRUE, dim_method = "pca", cluster_method = "leiden", show_plot = FALSE, return_plot = FALSE, save_plot = FALSE, ...) {
  if (!class(obj)[[1]] == "giotto") {
    stop("Object needs to be a Giotto object!")
  }

  if (!dim_method %in% c("pca", "tsne", "umap")) {
    message("Starting Giotto Workflow")
  }

  obj <- normalizeGiotto(obj)

  # HVF
  if (calculateHVG) {
    message("Giotto: Calculating HVFs")
    obj <- Giotto::calculateHVF(obj, show_plot = show_plot, return_plot = return_plot, save_plot = save_plot)
  }


  # Dimension Reduction
  message("Giotto: Performing Dimension Reduction")
  obj <- Giotto::runPCA(obj)

  obj <- switch(dim_method,
    pca = {
      obj
    },
    tsne = {
      Giotto::runtSNE(obj)
    },
    umap = {
      Giotto::runUMAP(obj)
    }
  )

  # Create Network using the performed dimension reduction
  message("Giotto: Calculating Network")
  obj <- Giotto::createNearestNetwork(
    gobject = obj,
    dim_reduction_to_use = dim_method,
    dim_reduction_name = dim_method
  )

  # Clustering
  message("Giotto: Clustering")
  obj <- switch(cluster_method,
    leiden = {
      Giotto::doLeidenCluster(obj)
    },
    louvain = {
      Giotto::doLouvainCluster(obj)
    },
    kmeans = {
      Giotto::doKmeans(obj,
        dim_reduction_to_use = dim_method,
        dim_reduction_name = dim_method
      )
    },
    hclust = {
      Giotto::doHclust(obj,
        dim_reduction_to_use = dim_method,
        dim_reduction_name = dim_method
      )
    },
  )

  message("Giotto: Finished")
  return(obj)
}
