#' Preprocess SingleCellExperiments and Spatial Experiments
#'
#' @param object SingleCellExperiment or SpatialExperiment, if AnnData or Seurat it will be converted
#' @param min_umi minimum umi count for spots/cells
#' @param assay assay to use for calculation, you can use any assay but counts is recommended
#'
#' @export
preprocess <- function(object, min_umi = 500, assay = "counts") {
  if (is.null(object)) {
    stop("Please provide an object")
  }

  # when not spatial convert to SCE
  if (!is(object, "SpatialExperiment")) {
    object <- convert_to_sce(object)
  }


  # Filtering
  # min UMI Count per Observation
  message("Removing ", sum(colSums(SummarizedExperiment::assay(object, assay)) < min_umi), " observations with umi count below threshold")
  object <- object[, colSums(SummarizedExperiment::assay(object, assay)) >= min_umi]

  # remove all zero genes
  message("Removing ", sum(rowSums(SummarizedExperiment::assay(object, assay)) == 0), " variables with all zero expression")
  object <- object[rowSums(SummarizedExperiment::assay(object, assay)) > 0, ]

  return(object)
}
