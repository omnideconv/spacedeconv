#' No signature calculated, just call the deconvolute method
#' @returns NULL

build_model_card <- function() {
  message("This method does not build a signature, just call the deconvolute method")
  return(NULL)
}

#' Deconvolute CARD
#'
#' @param sce Single Cell Experiment
#' @param spe Spatial Experiment
#' @param cell_type_col column of SCE containing cell type information
#' @param assay_sc assay of sce to use
#' @param assay_sp assay of spe to use
#' @param batch_id_col batch id column in spatialExperiment
deconvolute_card <- function(sce, spe, cell_type_col = "cell_ontology_class", assay_sc = "counts", assay_sp = "counts", batch_id_col = NULL) {


  # checks
  if (is.null(sce)) {
    stop("Parameter 'sce' is missing or null, but is required!")
  }

  if (is.null(spe)) {
    stop("Parameter 'spe' is missing or null, but is required")
  }

  # spatial expression
  spExpression <- SummarizedExperiment::assay(spe, assay_sp) %>% as("dgCMatrix")

  # spatial coords
  spCoords <- SpatialExperiment::spatialCoords(spe) %>% as.data.frame()
  colnames(spCoords) <- c("x", "y")

  # single cell expression
  scExpression <- SummarizedExperiment::assay(sce, assay_sc) %>% as("dgCMatrix")

  # create Metadata
  if (!is.null(batch_id_col)) {
    batchIDs <- SingleCellExperiment::colData(sce)[[batch_id_col]]
    scMeta <- data.frame(
      cellType = SingleCellExperiment::colData(sce)[[cell_type_col]],
      batchIDs = SingleCellExperiment::colData(sce)[[batch_id_col]],
      row.names = colnames(sce)
    )
    sample.varname <- "batchIDs"
  } else {
    scMeta <- data.frame(
      cellType = SingleCellExperiment::colData(sce)[[cell_type_col]],
      row.names = colnames(sce)
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


  return(deconvolution)
}
