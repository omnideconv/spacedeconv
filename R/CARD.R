#' No signature calculated, just call the deconvolute method
#' @returns NULL

build_model_card <- function() {
  message("This method does not build a signature, just call the deconvolute method")
  return(NULL)
}

#' Deconvolute CARD
#'
#' @param single_cell_obj Single Cell Experiment
#' @param spatial_obj Spatial Experiment
#' @param cell_type_col column of SCE containing cell type information
#' @param assay_sc assay of single_cell_obj to use
#' @param assay_sp assay of spatial_obj to use
#' @param batch_id_col batch id column in spatialExperiment
#' @param result_name token to identify deconvolution results in object, default = "card"
deconvolute_card <- function(single_cell_obj, spatial_obj, cell_type_col = "cell_ontology_class", assay_sc = "counts", assay_sp = "counts", batch_id_col = "orig.ident", result_name = "card") {
  # checks
  if (is.null(single_cell_obj)) {
    stop("Parameter 'single_cell_obj' is missing or null, but is required!")
  }

  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is missing or null, but is required")
  }

  if (is.null(batch_id_col) || !batch_id_col %in% names(colData(spatial_obj))){
    stop("Paramter 'batch_id_col' is missing, null or not available in single cell object")
  }

  # check if requested assay exists
  if (!assay_sc %in% names(SummarizedExperiment::assays(single_cell_obj))) {
    message(
      "requested assay ",
      assay_sc,
      " not available in single cell object. Using first available assay."
    )
    assay_sc <- names(SummarizedExperiment::assays(single_cell_obj))[1] # change to first available assay request not available
  }

  if (!assay_sp %in% names(SummarizedExperiment::assays(spatial_obj))) {
    message(
      "requested assay ",
      assay_sp,
      " not available in spatial object. Using first available assay."
    )
    assay_sp <- names(SummarizedExperiment::assays(spatial_obj))[1]
  }

  # spatial expression
  spExpression <- SummarizedExperiment::assay(spatial_obj, assay_sp) %>% as("dgCMatrix")

  # spatial coords
  spCoords <- SpatialExperiment::spatialCoords(spatial_obj) %>% as.data.frame()
  colnames(spCoords) <- c("x", "y")

  # single cell expression
  scExpression <- SummarizedExperiment::assay(single_cell_obj, assay_sc) %>% as("dgCMatrix")

  # create Metadata
  if (!is.null(batch_id_col)) {
    batchIDs <- SingleCellExperiment::colData(single_cell_obj)[[batch_id_col]]
    scMeta <- data.frame(
      cellType = SingleCellExperiment::colData(single_cell_obj)[[cell_type_col]],
      batchIDs = SingleCellExperiment::colData(single_cell_obj)[[batch_id_col]],
      row.names = colnames(single_cell_obj)
    )
    sample.varname <- "batchIDs"
  } else {
    scMeta <- data.frame(
      cellType = SingleCellExperiment::colData(single_cell_obj)[[cell_type_col]],
      row.names = colnames(single_cell_obj)
    )
    sample.varname <- NULL
  }

  # create CARD object
  cardObj <- CARD::createCARDObject(
    sc_count = scExpression,
    sc_meta = scMeta,
    spatial_count = spExpression,
    spatial_location = spCoords,
    sample.varname = sample.varname,
    ct.varname = "cellType",
    ct.select = NULL, # sure?
    minCountGene = 0, # sure?
    minCountSpot = 0 # sure?
  )

  # deconvolute
  cardObj <- CARD::CARD_deconvolution(cardObj)

  # extract results
  deconvolution <- cardObj@Proportion_CARD
  colnames(deconvolution) <- make.names(colnames(deconvolution))

  # attach token
  deconvolution <- attachToken(deconvolution, result_name)


  return(deconvolution)
}
