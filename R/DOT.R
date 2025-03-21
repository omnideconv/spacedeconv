#' No signature calculated, just call the deconvolute method
#' @returns NULL
build_model_dot <- function() {
  message("This method does not build a signature, just call the deconvolute method")
  return(NULL)
}



#' Deconvolute DOT
#'
#' @param single_cell_obj Single Cell Experiment
#' @param spatial_obj Spatial Experiment
#' @param cell_type_col column of SCE containing cell type information
#' @param assay_sc assay of single_cell_obj to use
#' @param assay_sp assay of spatial_obj to use
#' @param result_name token to identify deconvolution results in object, default = "dot"
#' @param ratios_weight This parameter determines how much the cell type ratios in the reference sc dataset are expected to match the ratios in the target spatial data. A higher value would match the cell type ratios more closely. By default this is 0 since we do not assume that the modalities match, but usually a value between 0.1 and 0.25 is reasonable depending on the quality of the sc data. If the two modalities are expected to match perfectly (e.g., come from the same tissue), a value closer to 1 would give the best results.
#' @param max_spot_size this parameter determines the maximum number of cells per spot and a higher value indicates lower resolution in the spatial data. The results should be robust for the default value but will give better estimates of the absolute abundances if the appropriate parameter is passed. For Visium, this is set to 20 (default). For ST it should be set to a higher value (e.g., 200).
#' @param ... additional parameters passed to DOTr methods
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
