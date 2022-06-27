#' List of supported deconvolution methods
#' @export
#'
deconvolution_methods <- c("RCTD" = "rctd", "DWLS" = "dwls")




#' Deconvolution
#' @param spatial_object A SpatialExperiment
#' @param single_cell_object A SingleCellExperiment
#' @param cell_type_col Column name of the single_cell_object where the cell type can be found
#' @param method Deconvolution Method to use, see deconvolution_methods() for a full list of available methods
#' @param ... Further parameters passed to the selected deconvolution method
#' @returns The deconvolution result as a table
#' @export
deconvolute <- function(spatial_object, single_cell_object, cell_type_col = "cell_ontology_class", method = NULL, ...) {
  if (is.null(spatial_object)) {
    stop("Parameter 'spatial_object' is missing or null, but is required.")
  }
  if (is.null(single_cell_object)) {
    stop("Parameter 'single_cell_object' is missing or null, but is required.")
  }

  # if got the methods name and not the token
  if (method %in% names(deconvolution_methods)) {
    method <- deconvolution_methods[[method]]
  }
  method <- tolower(method)

  # TODO Type checks for the spatial and single cell object
  # General type checks will be performed here, also matrix + annotation handling
  # Method specific processing steps will be located in the switch statement

  deconv <- switch(method,
    rctd = {
      deconvolute_rctd(...)
    },
    dwls = {
      deconvolute_omnideconv(method="dwls", ...)
    }
  )
  return(deconv)
}
