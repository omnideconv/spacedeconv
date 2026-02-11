#' No Signature for DOT
#'
#' DOT builds its model internally, so this function returns `NULL` and you
#' should call `deconvolute()` with the DOT method (see
#' `spacedeconv::deconvolution_methods`).
#'
#' @return NULL
build_model_dot <- function() {
  message("This method does not build a signature, just call the deconvolute method")
  return(NULL)
}



#' Deconvolute with DOT
#'
#' Runs DOT deconvolution and returns the DOT weight matrix with column names
#' prefixed by `result_name`.
#'
#' @param single_cell_obj `SingleCellExperiment`.
#' @param spatial_obj `SpatialExperiment`.
#' @param cell_type_col Column in `single_cell_obj` containing cell type labels.
#' @param assay_sc Single-cell assay to use.
#' @param assay_sp Spatial assay to use.
#' @param result_name Prefix used to label result columns (default: "dot").
#' @param ratios_weight Strength of matching cell-type ratios between scRNA-seq
#' and spatial data (default: 0).
#' @param max_spot_size Maximum number of cells per spot (platform dependent).
#' @param ... Additional parameters passed to DOT setup/run functions.
deconvolute_dot <- function(single_cell_obj,
                            spatial_obj,
                            cell_type_col = "cell_ontology_class",
                            assay_sc = "counts", assay_sp = "counts",
                            result_name = "dot",
                            ratios_weight = 0,
                            max_spot_size = 20,
                            ...) {
  if (is.null(single_cell_obj)) {
    stop("Parameter 'single_cell_obj' is missing or null, but is required!")
  }

  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is missing or null, but is required")
  }

  # check if requested SC assay exists
  if (!assay_sc %in% names(SummarizedExperiment::assays(single_cell_obj))) {
    message(
      "requested assay ",
      assay_sc,
      " not available in single cell object. Using first available assay."
    )
    assay_sc <- names(SummarizedExperiment::assays(single_cell_obj))[1] # change to first available assay request not available
  }

  # check if requested SP assay exists
  if (!assay_sp %in% names(SummarizedExperiment::assays(spatial_obj))) {
    message(
      "requested assay ",
      assay_sp,
      " not available in spatial object. Using first available assay."
    )
    assay_sp <- names(SummarizedExperiment::assays(spatial_obj))[1]
  }

  # spatial expression
  spExpression <- SummarizedExperiment::assay(
    spatial_obj, assay_sp
  ) %>% as("dgCMatrix")

  # spatial coords
  spCoords <- SpatialExperiment::spatialCoords(spatial_obj) %>% as.data.frame()
  colnames(spCoords) <- c("x", "y")

  # single cell expression and cell type labels
  scExpression <- SummarizedExperiment::assay(single_cell_obj, assay_sc) %>% as("dgCMatrix")
  cellTypes <- SingleCellExperiment::colData(single_cell_obj)[[cell_type_col]]

  # create DOT objects
  dot.srt <- setup.srt(srt_data = spExpression, srt_coords = spCoords, ...) # additional parameters passed as ...
  dot.ref <- setup.ref(ref_data = scExpression, ref_annotations = cellTypes, ...)

  dot <- create.DOT(dot.srt, dot.ref, ...)

  # Run Deconvolution
  dot <- run.DOT.lowresolution(dot, ratios_weight = ratios_weight, max_spot_size = max_spot_size, ...) # additional parameters passed as ...

  # attach token
  deconvolution_result <- attachToken(dot@weights, result_name) # add dot_ to column names

  return(deconvolution_result)
}
