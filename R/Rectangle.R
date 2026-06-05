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
#' @param assay_sp Spatial assay to use. Rectangle expects TPM-like values.
#' @param layer Optional AnnData layer for Rectangle.
#' @param raw Use AnnData raw values in Rectangle.
#' @param optimize_cutoffs Optimize Rectangle signature cutoffs.
#' @param p P-value cutoff used when `optimize_cutoffs = FALSE`.
#' @param lfc Log-fold-change cutoff used when `optimize_cutoffs = FALSE`.
#' @param n_cpus Number of CPU cores to use in Rectangle.
#' @param gene_expression_threshold Rectangle gene expression threshold.
#'
#' @return RectangleSignatureResult Python object.
build_model_rectangle <- function(single_cell_obj,
                                  spatial_obj = NULL,
                                  cell_type_col = "cell_ontology_class",
                                  assay_sc = "counts",
                                  assay_sp = "counts",
                                  layer = NULL,
                                  raw = FALSE,
                                  optimize_cutoffs = TRUE,
                                  p = 0.015,
                                  lfc = 1.5,
                                  n_cpus = NULL,
                                  gene_expression_threshold = 0.5) {
  if (is.null(single_cell_obj)) {
    stop("Parameter 'single_cell_obj' is missing or null, but is required.")
  }

  if (!checkCol(single_cell_obj, cell_type_col)) {
    stop(paste0("Column \"", cell_type_col, "\" can't be found in single cell object"))
  }

  ad <- rectangle_prepare_single_cell(single_cell_obj, cell_type_col, assay_sc)
  bulks <- NULL

  if (!is.null(spatial_obj)) {
    bulks <- rectangle_prepare_bulks(spatial_obj, assay_sp)
  }

  reticulate::source_python(system.file("python", "rectangle_spacedeconv.py", package = "spacedeconv"))

  signature <- py_build_rectangle_signatures(
    sc_obj = ad,
    bulks = bulks,
    cell_type_col = cell_type_col,
    layer = layer,
    raw = raw,
    optimize_cutoffs = optimize_cutoffs,
    p = p,
    lfc = lfc,
    n_cpus = n_cpus,
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
#' @param assay_sp Spatial assay to use. Rectangle expects TPM-like values.
#' @param result_name Prefix used to label result columns (default: "rectangle").
#' @param n_cpus Number of CPU cores to use in Rectangle.
#' @param correct_mrna_bias Correct mRNA bias in Rectangle.
#' @param ... Ignored. Present for compatibility with the generic dispatcher.
deconvolute_rectangle <- function(spatial_obj,
                                  signature = NULL,
                                  assay_sp = "counts",
                                  result_name = "rectangle",
                                  n_cpus = NULL,
                                  correct_mrna_bias = TRUE,
                                  ...) {
  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is missing or null, but is required.")
  }

  if (is.null(signature)) {
    stop("Parameter 'signature' is missing or null, but is required for Rectangle. Build it first with build_model().")
  }

  bulk_df <- rectangle_prepare_bulks(spatial_obj, assay_sp)

  reticulate::source_python(system.file("python", "rectangle_spacedeconv.py", package = "spacedeconv"))

  deconv <- py_deconvolute_rectangle(
    signature_result = signature,
    bulks = bulk_df,
    correct_mrna_bias = correct_mrna_bias,
    n_cpus = n_cpus
  )

  deconv <- as.matrix(deconv)
  if (nrow(deconv) == ncol(spatial_obj)) {
    rownames(deconv) <- colnames(spatial_obj)
  }
  deconv <- attachToken(deconv, result_name)

  return(deconv)
}

rectangle_prepare_single_cell <- function(single_cell_obj,
                                          cell_type_col,
                                          assay_sc) {
  if (!checkCol(single_cell_obj, cell_type_col)) {
    stop(paste0("Column \"", cell_type_col, "\" can't be found in single cell object"))
  }

  if (!assay_sc %in% names(SummarizedExperiment::assays(single_cell_obj))) {
    message(
      "requested assay ",
      assay_sc,
      " not available in single cell object. Using first available assay."
    )
    assay_sc <- names(SummarizedExperiment::assays(single_cell_obj))[1]
  }

  ad <- spe_to_ad(single_cell_obj, assay = assay_sc)
  ad$obs[[cell_type_col]] <- as.character(SingleCellExperiment::colData(single_cell_obj)[[cell_type_col]])
  ad$obs_names <- colnames(single_cell_obj)
  ad$var_names <- rownames(single_cell_obj)

  return(ad)
}

rectangle_prepare_bulks <- function(spatial_obj,
                                    assay_sp) {
  if (!assay_sp %in% names(SummarizedExperiment::assays(spatial_obj))) {
    message(
      "requested assay ",
      assay_sp,
      " not available in spatial object. Using first available assay."
    )
    assay_sp <- names(SummarizedExperiment::assays(spatial_obj))[1]
  }

  bulk_matrix <- Matrix::t(SummarizedExperiment::assay(spatial_obj, assay_sp))
  bulk_df <- as.data.frame(as.matrix(bulk_matrix))
  rownames(bulk_df) <- colnames(spatial_obj)
  colnames(bulk_df) <- rownames(spatial_obj)

  return(bulk_df)
}
