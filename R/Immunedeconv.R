#' No signature calculated, just call the deconvolute method
#' @returns NULL

build_model_immunedeconv <- function() {
  message("This method does not build a signature, just call the deconvolute method")
  return(NULL)
}


# is this located better somewhere else?
if (requireNamespace("xCell")){
  xCell.data <- xCell::xCell.data
}

#' Deconvolute Immundeconv
#' @param spe SpatialExperiment
#' @param method deconvolution algorithm
#' @param assay_sp assay of spatial object to use
#' @param ... further parameters passed to the selected method
deconvolute_immunedeconv <- function(spe, method = NULL, assay_sp = "counts", ...) {
  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if (is.null(method)) {
    stop("Parameter 'method' is missing or null, but is required")
  }

  # check for HGNC Symbols!!!
  rownames(spe) <- toupper(rownames(spe)) # temporary fix

  # extract bulk
  bulk <- NULL
  if (!is.null(spe)){
    bulk <- SummarizedExperiment::assay(spe, assay_sp)
    bulk <- as.matrix(bulk)
  }

  # some parameters are handled with the ... option
  deconv <- immunedeconv::deconvolute(
    gene_expression = bulk,
    method = method,
    indications = NULL,
    tumor = TRUE,
    arrays = FALSE
  )

  return(convertImmunedeconvMatrix(deconv))
}


#' Deconvolute Immunedeconv mouse
#' @param spe SpatialExperiment
#' @param method deconvolution algorithm
#' @param rmgenes genes to remove from the analysis
#' @param algorithm statistical algorithm for SeqImmuCC (ignored by all other methods)
#' @param assay_sp assay of spatial object to use
#' @param ... additional parameters passed to function
#'
deconvolute_immunedeconv_mouse <- function(spe, method = NULL, rmgenes = NULL, algorithm = NULL, assay_sp = "counts", ...) {
  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if (is.null(method)) {
    stop("Parameter 'method' is missing or null, but is required")
  }

  # extract bulk
  bulk <- NULL
  if (!is.null(spe)){
    bulk <- SummarizedExperiment::assay(spe, assay_sp)
    bulk <- as.matrix(bulk)
  }

  deconv <- immunedeconv::deconvolute_mouse(
    gene_expression_matrix = bulk,
    method = method,
    rmgenes = rmgenes,
    algorithm = algorithm
  )

  return(convertImmunedeconvMatrix(deconv))
}

#' Convert Immunedeconv Matrix to match SpaceDeconv Format
#' @param deconvResult immunedeconv result matrix
#' @returns transformed matrix where cell types are columns
convertImmunedeconvMatrix <- function(deconvResult) {
  # TODO Checks


  # use column 1 as rownames and transpose
  result <- deconvResult[, 2:ncol(deconvResult)]
  result <- as.matrix(result)
  rownames(result) <- deconvResult$cell_type
  result <- t(result)

  return(result)
}
