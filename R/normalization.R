#' Normalize Expression Object
#'
#' @param object SingleCellExperiment
#' @param method normalization method, ("cpm")
#' @export
normalize <- function(object, method = "cpm") {
  if (is.null(object)) {
    stop("Parameter 'object' is null or missing, but is required!")
  }

  # check if rownames and colnames are set
  if (checkRowColumn(single_cell_obj)||checkRowColumn(spatial_obj)){
    stop ("Rownames or colnames not set for single_cell_obj or spatial_obj but need to be available!")
  }

  # ensure library size > 0
  object <- removeZeroExpression(object)

  if (class(object)[[1]] %in% c("SingleCellExperiment", "SpatialExperiment")) {
    if (method == "cpm") {
      SummarizedExperiment::assay(object, "cpm") <- as(edgeR::cpm(object), "dgCMatrix")
    }
  } else {
    message("normalization currently only implemented for SingleCellExperiment and SpatialExperiment")
  }

  message("Normalized object using ", method, ". Note that the normalization is saved in an additional assay.")

  return(object)
}
