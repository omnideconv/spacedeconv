#' Preprocess SingleCellExperiments and Spatial Experiments
#'
#' @param object SingleCellExperiment or SpatialExperiment, if AnnData or Seurat it will be converted
#' @param min_umi minimum umi count for spots/cells
#' @param max_umi maximimum umi
#' @param assay assay to use for calculation, you can use any assay but counts is recommended
#'
#' @export
preprocess <- function(object, min_umi = 500, max_umi = NULL, assay = "counts") {
  cli::cli_rule(left = "spacedeconv")
  cli::cli_progress_step("testing parameter", msg_done = "parameter OK")

  if (is.null(object)) {
    stop("Please provide an object")
  }

  # when not spatial convert to SCE
  if (!is(object, "SpatialExperiment")) {
    object <- convert_to_sce(object)
  }

  # Filtering
  # min UMI Count per Observation
  nObservation <- sum(colSums(SummarizedExperiment::assay(object, assay)) < min_umi)
  cli::cli_progress_step(
    msg = paste0("Removing ", nObservation, " observations with umi count below threshold"),
    msg_done = paste0("Removed ", nObservation, " observations with umi count below threshold")
  )

  object <- object[, colSums(SummarizedExperiment::assay(object, assay)) >= min_umi]


  # max UMI Count per Observation, only if value is provided
  if (!is.null(max_umi)) {
    nObservation <- sum(colSums(SummarizedExperiment::assay(object, assay)) > max_umi)
    cli::cli_progress_step(
      msg = paste0("Removing ", nObservation, " observations with umi count abvove threshold"),
      msg_done = paste0("Removed ", nObservation, " observations with umi count above threshold")
    )

    object <- object[, colSums(SummarizedExperiment::assay(object, assay)) <= max_umi]
  }

  # remove all zero genes
  nVariable <- sum(rowSums(SummarizedExperiment::assay(object, assay)) == 0)
  cli::cli_progress_step(
    msg = paste0("Removing ", nVariable, " variables with all zero expression"),
    msg_done = paste0("Removed ", nVariable, " variables with all zero expression")
  )

  object <- object[rowSums(SummarizedExperiment::assay(object, assay)) > 0, ]

  cli::cli_process_done()

  return(object)
}
