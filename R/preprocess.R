#' Preprocess Single Cell and Spatial Data for spacedeconv
#'
#' Converts input to `SingleCellExperiment` when needed and applies common QC
#' filters: UMI thresholds, removal of all-zero genes, optional mitochondrial
#' gene removal, and duplicate-gene handling. Also emits a warning if Ensembl
#' identifiers are detected.
#'
#' @param object Input object (`SingleCellExperiment`, `SpatialExperiment`, AnnData,
#' or Seurat). AnnData/Seurat are converted to `SingleCellExperiment`.
#' @param min_umi Minimum UMI count per cell/spot to keep.
#' @param max_umi Maximum UMI count per cell/spot to keep (optional).
#' @param assay Assay to use for filtering (default: "counts").
#' @param remove_mito Logical; drop genes with `MT-` prefix.
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

  # Identify and process duplicate row names
  duplicated_row_names <- rownames(SummarizedExperiment::assay(object, assay))[duplicated(rownames(SummarizedExperiment::assay(object, assay)))]

  if (length(duplicated_row_names) > 0) {
    cli::cli_alert_warning("Duplicated genes found: ", toString(duplicated_row_names))

    cli::cli_progress_step(
      msg = "Selecting gene with highest mean expression for duplicated genes",
      msg_done = "Selected gene with highest mean expression for duplicated genes"
    )

    # For each duplicated gene, keep the one with the highest mean expression
    unique_genes <- lapply(duplicated_row_names, function(gene) {
      gene_rows <- which(rownames(SummarizedExperiment::assay(object, assay)) == gene)
      mean_expressions <- Matrix::rowMeans(SummarizedExperiment::assay(object, assay)[gene_rows, , drop = FALSE])
      gene_rows[which.max(mean_expressions)]
    })

    # Get all unique gene rows (non-duplicates + highest mean duplicates)
    all_gene_rows <- unique(c(which(!rownames(SummarizedExperiment::assay(object, assay)) %in% duplicated_row_names), unlist(unique_genes)))

    # Subset the object to keep only the selected rows
    object <- object[all_gene_rows, ]
  }



  cli::cli_progress_step(
    msg = "Checking for ENSEMBL Identifiers",
    msg_done = "Finished Preprocessing"
  )
  # Check for Ensembl identifiers in row names
  if (any(grepl("^ENS", rownames(SummarizedExperiment::assay(object, assay))))) {
    cli::cli_alert_warning("Warning: ENSEMBL identifiers detected in gene names")
    cli::cli_alert_info("Consider using Gene Names for first-generation deconvolution tools")
  }


  cli::cli_process_done()

  return(object)
}
