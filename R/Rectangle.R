#' Build Rectangle Signatures
#'
#' Builds a `RectangleSignatureResult` from single-cell count data. If a spatial
#' object is provided, Rectangle can use the matching bulk/spatial expression
#' profiles while optimizing signature cutoffs.
#'
#' @param single_cell_obj `SingleCellExperiment`.
#' @param spatial_obj Optional `SpatialExperiment`.
#' @param cell_type_col Column in `single_cell_obj` containing cell type labels.
#' @param assay_sc Single-cell assay to use. Rectangle expects raw counts.
#' @param assay_sp Spatial assay to use. Use raw UMI counts or CPM for UMI-based spatial data;
#' Rectangle rescales each spot to a total of one million. For non-UMI
#' full-length data, supply TPM. Do not use log-transformed values.
#' @param optimize_cutoffs Optimize Rectangle signature cutoffs.
#' @param p P-value cutoff used when `optimize_cutoffs = FALSE`.
#' @param lfc Log-fold-change cutoff used when `optimize_cutoffs = FALSE`.
#' @param n_cpus Number of CPU cores to use in Rectangle (default: 1).
#' @param gene_expression_threshold Rectangle gene expression threshold.
#'
#' @return RectangleSignatureResult Python object.
build_model_rectangle <- function(single_cell_obj,
                                  spatial_obj = NULL,
                                  cell_type_col = "cell_ontology_class",
                                  assay_sc = "counts",
                                  assay_sp = "counts",
                                  optimize_cutoffs = TRUE,
                                  p = 0.015,
                                  lfc = 1.5,
                                  n_cpus = 1L,
                                  gene_expression_threshold = 0.5) {
  if (is.null(single_cell_obj)) {
    stop("Parameter 'single_cell_obj' is missing or null, but is required.")
  }

  if (!checkCol(single_cell_obj, cell_type_col)) {
    stop(paste0("Column \"", cell_type_col, "\" can't be found in single cell object"))
  }

  rectangle_check_cpus(n_cpus)
  ad <- rectangle_prepare_single_cell(single_cell_obj, cell_type_col, assay_sc)
  bulks <- NULL

  if (!is.null(spatial_obj)) {
    bulks <- rectangle_prepare_bulks(spatial_obj, assay_sp)
  }

  python <- rectangle_python()

  signature <- python$py_build_rectangle_signatures(
    sc_obj = ad,
    bulks = bulks,
    cell_type_col = cell_type_col,
    optimize_cutoffs = optimize_cutoffs,
    p = p,
    lfc = lfc,
    n_cpus = as.integer(n_cpus),
    gene_expression_threshold = gene_expression_threshold
  )

  return(signature)
}

#' Deconvolute with Rectangle
#'
#' Runs Rectangle deconvolution via Python/reticulate and returns a matrix with
#' column names prefixed by `result_name`.
#'
#' @param spatial_obj `SpatialExperiment`.
#' @param signature RectangleSignatureResult from `build_model()`.
#' @param assay_sp Spatial assay to use. Use raw UMI counts or CPM for UMI-based spatial data;
#' Rectangle rescales each spot to a total of one million. For non-UMI
#' full-length data, supply TPM. Do not use log-transformed values.
#' @param result_name Prefix used to label result columns (default: "rectangle").
#' @param n_cpus Number of CPU cores to use in Rectangle (default: 1).
#' @param correct_mrna_bias Correct mRNA bias in Rectangle.
#' @param ... Ignored. Present for compatibility with the generic dispatcher.
deconvolute_rectangle <- function(spatial_obj,
                                  signature = NULL,
                                  assay_sp = "counts",
                                  result_name = "rectangle",
                                  n_cpus = 1L,
                                  correct_mrna_bias = TRUE,
                                  ...) {
  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is missing or null, but is required.")
  }

  if (is.null(signature)) {
    stop("Parameter 'signature' is missing or null, but is required for Rectangle. Build it first with build_model().")
  }

  rectangle_check_cpus(n_cpus)
  bulk_df <- rectangle_prepare_bulks(spatial_obj, assay_sp)

  python <- rectangle_python()

  deconv <- python$py_deconvolute_rectangle(
    signature_result = signature,
    bulks = bulk_df,
    correct_mrna_bias = correct_mrna_bias,
    n_cpus = as.integer(n_cpus)
  )

  deconv <- as.matrix(deconv)
  if (!is.numeric(deconv) || !all(is.finite(deconv)) || any(deconv < -1e-8)) {
    stop("Rectangle returned invalid cell-type estimates.")
  }
  if (is.null(rownames(deconv)) || anyDuplicated(rownames(deconv)) ||
      !setequal(rownames(deconv), colnames(spatial_obj))) {
    stop("Rectangle result spot IDs do not match the spatial object.")
  }
  deconv <- deconv[colnames(spatial_obj), , drop = FALSE]
  deconv <- attachToken(deconv, result_name)

  return(deconv)
}

rectangle_prepare_single_cell <- function(single_cell_obj,
                                          cell_type_col,
                                          assay_sc) {
  if (!checkCol(single_cell_obj, cell_type_col)) {
    stop(paste0("Column \"", cell_type_col, "\" can't be found in single cell object"))
  }

  rectangle_check_assay(single_cell_obj, assay_sc, counts = TRUE)
  labels <- SummarizedExperiment::colData(single_cell_obj)[[cell_type_col]]
  if (anyNA(labels) || any(!nzchar(as.character(labels))) || length(unique(labels)) < 2L) {
    stop("Rectangle requires non-missing cell type labels and at least two cell types.")
  }

  ad <- spe_to_ad(single_cell_obj, assay = assay_sc)
  ad$obs[[cell_type_col]] <- as.character(SingleCellExperiment::colData(single_cell_obj)[[cell_type_col]])
  ad$obs_names <- colnames(single_cell_obj)
  ad$var_names <- rownames(single_cell_obj)

  return(ad)
}

rectangle_prepare_bulks <- function(spatial_obj,
                                    assay_sp) {
  rectangle_check_assay(spatial_obj, assay_sp)

  bulk_matrix <- Matrix::t(SummarizedExperiment::assay(spatial_obj, assay_sp))
  bulk_df <- as.data.frame(as.matrix(bulk_matrix))
  rownames(bulk_df) <- colnames(spatial_obj)
  colnames(bulk_df) <- rownames(spatial_obj)

  return(bulk_df)
}

rectangle_check_assay <- function(object, assay, counts = FALSE) {
  if (length(assay) != 1L || is.na(assay) || !assay %in% SummarizedExperiment::assayNames(object)) {
    stop("Requested Rectangle assay is not available: ", assay)
  }
  for (ids in list(rownames(object), colnames(object))) {
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
      stop("Rectangle requires unique, non-empty gene and cell/spot IDs.")
    }
  }
  x <- SummarizedExperiment::assay(object, assay)
  if (!all(is.finite(x)) || any(x < 0)) {
    stop("Rectangle requires finite, non-negative expression values.")
  }
  if (counts && any(x != round(x))) {
    stop("Rectangle single-cell reference requires raw integer counts.")
  }
  if (any(colSums(x) <= 0)) {
    stop("Rectangle requires positive library sizes for every cell/spot.")
  }
}

rectangle_check_cpus <- function(n_cpus) {
  if (!is.numeric(n_cpus) || length(n_cpus) != 1L || !is.finite(n_cpus) ||
      n_cpus < 1 || n_cpus != round(n_cpus)) {
    stop("n_cpus must be a positive integer.")
  }
}

rectangle_python <- function() {
  if (!reticulate::py_module_available("rectanglepy")) {
    stop("Rectangle requires rectanglepy in the active Python environment. See the installation instructions in README.md.")
  }
  python <- new.env(parent = baseenv())
  reticulate::source_python(system.file("python", "rectangle_spacedeconv.py",
    package = "spacedeconv", mustWork = TRUE), envir = python)
  python
}
