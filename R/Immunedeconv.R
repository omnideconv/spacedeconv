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
  deconv = as.data.frame(deconv)
  rownames(deconv) = deconv$cell_type
  deconv = deconv[, -1]
  deconv = t(deconv)
  deconv = as.matrix(deconv) # ???

  return(deconv)

}


#' Deconvolute Immunedeconv mouse
#' @param spe SpatialExperiment
#' @param method deconvolution algorithm
#' @param ... additional parameters passed to function
#'
deconvolute_immunedeconv_mouse <- function (spe, method=NULL, ...){
  if (is.null(spe)){
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if(is.null(method)){
    stop("Parameter 'method' is missing or null, but is required")
  }

  # extract bulk
  bulk <- SingleCellExperiment::counts(spe)
  bulk <- as.matrix(bulk)

  deconv = immunedeconv::deconvolute_mouse(bulk, method = method)
}
