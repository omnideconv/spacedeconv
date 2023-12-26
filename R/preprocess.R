#' Preprocess SingleCellExperiments and Spatial Experiments
#'
#' @param object SingleCellExperiment or SpatialExperiment, if AnnData or Seurat it will be converted
#' @param min_umi minimum umi count for spots/cells
#' @param max_umi maximimum umi
#' @param assay assay to use for calculation, you can use any assay but counts is recommended
#' @param remove_mito remove mitochondria genes
#'
#' @export
preprocess <- function(object, min_umi = 500, max_umi = NULL, assay = "counts", remove_mito = FALSE) {
  cli::cli_rule(left = "spacedeconv")
  cli::cli_progress_step("testing parameter", msg_done = "parameter OK")

  if (is.null(object)) {
    stop("Please provide an object")
  }

  # when not spatial convert to SCE
  if (!is(object, "SpatialExperiment")) {
    object <- convert_to_sce(object)
  }

  # convert to sparse matrices
  object <- check_datatype(object)

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

  # Check if mitochondrial genes removal is requested
  if (remove_mito) {
    mito_genes <- grep("^MT-", rownames(SummarizedExperiment::assay(object, assay)), value = TRUE)
    nMito <- length(mito_genes)

    # If mitochondrial genes are present, remove them
    if (nMito > 0) {
      cli::cli_progress_step(
        msg = paste0("Removing ", nMito, " mitochondria genes"),
        msg_done = paste0("Removed ", nMito, " mitochondria genes")
      )
      object <- object[!rownames(SummarizedExperiment::assay(object, assay)) %in% mito_genes, ]
    }
  } else {
    # Check for mitochondrial genes presence when removal is not requested
    mito_genes <- grep("^MT-", rownames(SummarizedExperiment::assay(object, assay)), value = TRUE)
    nMito <- length(mito_genes)

    # Display a warning if mitochondrial genes are present
    if (nMito > 0) {
      warning("There are ", nMito, " mitochondrial genes present. Consider removing them.")
    }
  }

  # Identify and print duplicate row names
  duplicated_row_names <- rownames(SummarizedExperiment::assay(object, assay))[duplicated(rownames(SummarizedExperiment::assay(object, assay)))]
  if (length(duplicated_row_names) > 0) {
    cli::cli_alert_warning("Duplicated genes found: ", toString(duplicated_row_names))
  }

  cli::cli_progress_step(
    msg = "Removing duplicated genes",
    msg_done = "Removed duplicated genes"
  )

  # Remove duplicate row names
  unique_row_names <- !duplicated(rownames(SummarizedExperiment::assay(object, assay)))
  object <- object[unique_row_names, ]


  cli::cli_progress_step(msg = "Checking for ENSEMBL Identifiers",
                         msg_done = "Finished Preprocessing")
  # Check for Ensembl identifiers in row names
  if (any(grepl("^ENS", rownames(SummarizedExperiment::assay(object, assay))))) {
    cli::cli_alert_warning("Warning: ENSEMBL identifiers detected in gene names")
    cli::cli_alert_info("Consider using Gene Names for first-generation deconvolution tools")
  }


  cli::cli_process_done()

  return(object)
}
