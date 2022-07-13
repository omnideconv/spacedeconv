#' Normalize Expression Object
#'
#' @param object SingleCellExperiment
#' @param method normalization method, ("cpm")
#' @export
normalize <- function(object, method = "cpm") {
  if (is.null(object)) {
    stop("Parameter 'object' is null or missing, but is required!")
  }

  if (class(object)[[1]] %in% c("SingleCellExperiment", "SpatialExperiment")) {
    if (method == "cpm") {
      SingleCellExperiment::cpm(object) <- edgeR::cpm(object)
    }
  } else {
    message("normalization currently only implemented for SingleCellExperiment and SpatialExperiment")
  }

  return(object)
}
