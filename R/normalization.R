#' Normalize Expression Data in SCE/SPE Objects
#'
#' Adds a normalized assay (`cpm`, `logcpm`, or `logpf`) to a
#' `SingleCellExperiment` or `SpatialExperiment`. The original assay
#' is left unchanged. `cpm` computes counts per million, `logcpm` is
#' `log(cpm + 1)`, and `logpf` is `log((counts + 1) / (gene_mean + 1))`
#' with gene-wise means computed from `counts + 1`.
#'
#' @param object A `SingleCellExperiment` or `SpatialExperiment` with raw counts.
#' @param method Normalization method: `"cpm"`, `"logcpm"`, or `"logpf"`.
#' @param assay Assay to normalize (default: "counts").
#'
#' @export
normalize <- function(object, method = "cpm", assay = "counts") {
  cli::cli_rule(left = "spacedeconv")

  cli::cli_progress_step("testing parameter", msg_done = "parameter OK")


  if (is.null(object)) {
    stop("Parameter 'object' is null or missing, but is required!")
  }

  # check if rownames and colnames are set
  if (checkRowColumn(object)) {
    stop("Rownames or colnames not set for expression object but need to be available!")
  }

  # convert to sparse matrices
  object <- check_datatype(object)

  cli::cli_progress_step(msg = paste0("Normalizing using ", method), msg_done = paste0("Finished normalization using ", method))

  if (class(object)[[1]] %in% c("SingleCellExperiment", "SpatialExperiment")) {
    if (method == "cpm") {
      SummarizedExperiment::assay(object, "cpm") <- as(edgeR::cpm(SummarizedExperiment::assay(object, assay)), "dgCMatrix")
    } else if (method == "logcpm") {
      SummarizedExperiment::assay(object, "logcpm") <- as(log(edgeR::cpm(SummarizedExperiment::assay(object, assay)) + 1), "dgCMatrix") # log(cpm+1)
    } else if (method == "logpf") {
      assay_data <- SummarizedExperiment::assay(object, assay) + 1
      mean_read_count <- Matrix::rowMeans(assay_data)
      SummarizedExperiment::assay(object, "logpf") <- as(log(assay_data / (mean_read_count + 1)), "dgCMatrix")
    }
  } else {
    message("normalization currently only implemented for SingleCellExperiment and SpatialExperiment")
  }

  cli::cli_progress_done()

  cli::cli_alert_info("Please note the normalization is stored in an additional assay")

  return(object)
}
