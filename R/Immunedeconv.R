#' No signature calculated, just call the deconvolute method
#' @returns NULL

build_model_immunedeconv <- function (){
  message("This method does not build a signature, just call the deconvolute method")
  return (NULL)
}


#' Deconvolute Immundeconv
#' @param spe SpatialExperiment
#' @param method deconvolution algorithm
#' @param ... further parameters passed to the selected method
deconvolute_immunedeconv <- function(spe, method = NULL, ...){

  if (is.null(spe)){
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if(is.null(method)){
    stop("Parameter 'method' is missing or null, but is required")
  }

  # check for HGNC Symbols!!!
  # this is only
  rownames(spe) = toupper(rownames(spe)) # temporary fix

  # extract bulk
  bulk <- SingleCellExperiment::counts(spe)
  bulk <- as.matrix(bulk)

  # some parameters are handled with the ... option
  deconv = immunedeconv::deconvolute(
    gene_expression = bulk,
    method=method,
    indications =NULL,
    tumor = TRUE,
    arrays=FALSE)

  # convert data to matrix where cell types are rows and samples are columns
  deconv = as.matrix(t(deconv)) # samples as rows
  rownames(deconv) <- deconv[, 1] # turn first column to rownames
  deconv = deconv[, -1] # remove first column containing cell types

  return(deconv)

}
