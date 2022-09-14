#' Normalize Expression Object
#'
#' @param object SingleCellExperiment
#' @param method normalization method, ("cpm")
#' @export
normalize <- function(object, method = "cpm") {
  if (is.null(object)) {
    stop("Parameter 'object' is null or missing, but is required!")
  }

  # ensure that library size > 0
  nspots = sum (Matrix::colSums(counts(object))==0)
  if (nspots>0){
    # remove spots with all zero expression
    message ("removing ", nspots, " spots with zero expression")
    object <- object[, !Matrix::colSums(counts(object))==0]
  }

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
