#' Normalize Expression Object
#'
#' @param object SingleCellExperiment
#' @param method normalization method, ("cpm", "logcpm")
#' @param assay which assay to use
#' @export
normalize <- function(object, method = "cpm", assay="counts") {
  cli::cli_rule(left = "spacedeconv")

  cli::cli_progress_step("testing parameter", msg_done = "parameter OK")


  if (is.null(object)) {
    stop("Parameter 'object' is null or missing, but is required!")
  }

  # check if rownames and colnames are set
  if (checkRowColumn(object)) {
    stop("Rownames or colnames not set for expression object but need to be available!")
  }

  cli::cli_progress_step(msg = paste0("Normalizing using ", method), msg_done = paste0("Finished normalization using ", method))

  if (class(object)[[1]] %in% c("SingleCellExperiment", "SpatialExperiment")) {
    if (method == "cpm") {
      SummarizedExperiment::assay(object, "cpm") <- as(edgeR::cpm(SummarizedExperiment::assay(object, assay)), "dgCMatrix")
    } else if (method == "logcpm") {
      SummarizedExperiment::assay(object, "logcpm") <- as(log(edgeR::cpm(SummarizedExperiment::assay(object, assay)) + 1), "dgCMatrix") # log(cpm+1)
    }
  } else {
    message("normalization currently only implemented for SingleCellExperiment and SpatialExperiment")
  }

  cli::cli_progress_done()

  cli::cli_alert_info("Please note the normalization is stored in an additional assay")

  return(object)
}
