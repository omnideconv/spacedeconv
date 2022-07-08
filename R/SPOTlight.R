#' Build Model SPOTlight
#'
#' @param sce SingleCellExperiment
#' @param cell_type_col Column where cell type information can be found
#' @param spe SpatialExperiment
#' @param markers (Optional) Marker Gene DataFrame, if NULL markers will be calculated from 'sce' based on the authors suggestion
build_model_spotlight <- function(sce, cell_type_col = "cell_ontology_class", spe, markers = NULL) {
  if (is.null(sce)) {
    stop("Parameter 'sce' is null or missing, but is required")
  }

  if (!checkCol(sce, cell_type_col)) {
    stop(paste0("Column \"", cell_type_col, "\" can't be found in single cell object"))
  }

  groups <- colData(sce)[[cell_type_col]] # cell type vector
  if (is.null(markers)) {
    message("No markers provided, calculating markers based on the authors suggestion")
    mgs <- getMarkersSPOTlight(
      sce = sce,
      cell_type_col = cell_type_col
    )
  }

  model <- SPOTlight::trainNMF(
    x = sce,
    y = spe,
    groups = groups,
    mgs = mgs,
    weight_id = "mean.AUC"
  )

  return(model)
}


#' Deconvolute SPOTlight
#'
#' @param spe SpatialExperiment
#' @param model SPOTlight Model
deconvolute_spotlight <- function(spe, model = NULL) {
  if (is.null(spe)) {
    stop("Parameter 'spe' missing or null, but is required")
  }

  if (is.null(model)) {
    stop("Model is missing or null, but is required")
  }

  # extract model information
  mod <- model$mod
  ref <- model$topic

  # deconvolute
  deconv <- SPOTlight::runDeconvolution(
    x = spe,
    mod = mod,
    ref = ref
  )
  return(deconv$mat)
}

#' Calculate Markers
#'
#' @param sce SingleCellExperiment
#' @param cell_type_col Column containing the cell type
#'
#' This Procedure reflects the suggestions of the SPOTlight authors, however,
#' they also state that there are other ways to calculate markers
getMarkersSPOTlight <- function(sce, cell_type_col = "cell_ontology_class") {
  # TODO checks

  groups <- colData(sce)[[cell_type_col]]
  sce <- scuttle::logNormCounts(sce) #  only if not log normalized yet!?
  mgs <- scran::scoreMarkers(sce, groups = groups)
  mgs_fil <- lapply(names(mgs), function(i) {
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.8, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
  })
  mgs_df <- do.call(rbind, mgs_fil)

  return(mgs_df) # could also return a subset
}
