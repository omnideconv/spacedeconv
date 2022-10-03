#' Build Model SPOTlight
#'
#' @param single_cell_obj SingleCellExperiment
#' @param cell_type_col Column where cell type information can be found
#' @param spatial_obj SpatialExperiment
#' @param assay_sc single cell assay to use
#' @param assay_sp spatial assay to use
#' @param markers (Optional) Marker Gene DataFrame, if NULL markers will be calculated from 'single_cell_obj' based on the authors suggestion
build_model_spotlight <- function(single_cell_obj, cell_type_col = "cell_ontology_class", spatial_obj, assay_sc = "counts", assay_sp = "counts", markers = NULL) {
  if (is.null(single_cell_obj)) {
    stop("Parameter 'single_cell_obj' is null or missing, but is required")
  }

  if (!checkCol(single_cell_obj, cell_type_col)) {
    stop(paste0("Column \"", cell_type_col, "\" can't be found in single cell object"))
  }

  if (is.null(spatial_obj)){
    stop("Parameter 'spatial'obj' is null or missing, but is required")
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

  # check if requested assay exists
  if (!assay_sp %in% names(SummarizedExperiment::assays(spatial_obj))) {
    message(
      "requested assay ",
      assay_sp,
      " not available in expression object. Using first available assay."
    )
    assay_sp <- names(SummarizedExperiment::assays(spatial_obj))[1] # change to first available assay request not available
  }

  groups <- colData(single_cell_obj)[[cell_type_col]] # cell type vector
  if (is.null(markers)) {
    message("No markers provided, calculating markers based on the authors suggestion")
    mgs <- getMarkersSPOTlight(
      single_cell_obj = single_cell_obj,
      cell_type_col = cell_type_col
    )
  }

  model <- SPOTlight::trainNMF(
    x = single_cell_obj,
    y = spatial_obj,
    groups = groups,
    mgs = mgs,
    weight_id = "mean.AUC",
    slot_sc = assay_sc, # not sure about this one!
    slot_sp = assay_sp
  )

  return(model)
}


#' Deconvolute SPOTlight
#'
#' @param spatial_obj SpatialExperiment
#' @param model SPOTlight Model
#' @param assay_sp spatial assay to use for the computation
#' @param result_name token to identify deconvolution results in object, default = "spotlight"
deconvolute_spotlight <- function(spatial_obj, model = NULL, assay_sp = "counts", result_name = "spotlight") {
  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' missing or null, but is required")
  }

  if (is.null(model)) {
    stop("Model is missing or null, but is required")
  }

  # check for model integrity: names(signature) must be "mod" "topic"

  # extract model information
  mod <- model$mod
  ref <- model$topic

  # deconvolute
  deconv <- SPOTlight::runDeconvolution(
    x = spatial_obj,
    mod = mod,
    ref = ref,
    slot = assay_sp
  )

  deconvolution <- deconv$mat

  # attach method token
  deconvolution <- attachToken(deconvolution, result_name)

  return(deconvolution)
}

#' Calculate Markers
#'
#' @param single_cell_obj SingleCellExperiment
#' @param cell_type_col Column containing the cell type
#'
#' This Procedure reflects the suggestions of the SPOTlight authors, however,
#' they also state that there are other ways to calculate markers
getMarkersSPOTlight <- function(single_cell_obj, cell_type_col = "cell_ontology_class") {
  # TODO checks! Check if cell_type_col actually exists in single_cell_obj

  groups <- colData(single_cell_obj)[[cell_type_col]]
  single_cell_obj <- scuttle::logNormCounts(single_cell_obj) #  only if not log normalized yet!?
  mgs <- scran::scoreMarkers(single_cell_obj, groups = groups)
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
