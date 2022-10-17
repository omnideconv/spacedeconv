#' Normalize Expression Object
#'
#' @param object SingleCellExperiment
#' @param method normalization method, ("cpm", "logcpm")
#' @export
normalize <- function(object, method = "cpm") {
  if (is.null(object)) {
    stop("Parameter 'object' is null or missing, but is required!")
  }

  # check if rownames and colnames are set
  if (checkRowColumn(object)) {
    stop("Rownames or colnames not set for expression object but need to be available!")
  }

  # ensure library size > 0
  object <- removeZeroExpression(object)

  if (class(object)[[1]] %in% c("SingleCellExperiment", "SpatialExperiment")) {
    if (method == "cpm") {
      SummarizedExperiment::assay(object, "cpm") <- as(edgeR::cpm(object), "dgCMatrix")
    } else if (method == "logcpm") {
      SummarizedExperiment::assay(object, "logcpm") <- as(log(edgeR::cpm(object) + 1), "dgCMatrix") # log(cpm+1)
    }
  } else {
    message("normalization currently only implemented for SingleCellExperiment and SpatialExperiment")
  }

  message("Normalized object using ", method, ". Note that the normalization is saved in an additional assay.")

  return(object)
}
