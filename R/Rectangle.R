#' No Signature for Rectangle
#'
#' Rectangle builds its model internally, so this function returns `NULL` and
#' you should call `deconvolute()` with the Rectangle method.
#'
#' @return NULL
build_model_rectangle <- function() {
  message("Rectangle builds its model internally, please just use the deconvolute method")
  return(NULL)
}

#' Deconvolute with Rectangle
#'
#' Runs Rectangle deconvolution via Python/reticulate and returns a matrix with
#' column names prefixed by `result_name`.
#'
#' @param single_cell_obj `SingleCellExperiment`.
#' @param spatial_obj `SpatialExperiment`.
#' @param cell_type_col Column in `single_cell_obj` containing cell type labels.
#' @param assay_sc Single-cell assay to use. Rectangle expects raw counts.
#' @param assay_sp Spatial assay to use. Rectangle expects TPM-like values.
#' @param result_name Prefix used to label result columns (default: "rectangle").
#' @param n_cpus Number of CPU cores to use in Rectangle.
#' @param layer Optional AnnData layer for Rectangle.
#' @param raw Use AnnData raw values in Rectangle.
#' @param correct_mrna_bias Correct mRNA bias in Rectangle.
#' @param optimize_cutoffs Optimize Rectangle signature cutoffs.
#' @param p P-value cutoff used when `optimize_cutoffs = FALSE`.
#' @param lfc Log-fold-change cutoff used when `optimize_cutoffs = FALSE`.
#' @param gene_expression_threshold Rectangle gene expression threshold.
deconvolute_rectangle <- function(single_cell_obj,
                                  spatial_obj,
                                  cell_type_col = "cell_ontology_class",
                                  assay_sc = "counts",
                                  assay_sp = "counts",
                                  result_name = "rectangle",
                                  n_cpus = NULL,
                                  layer = NULL,
                                  raw = FALSE,
                                  correct_mrna_bias = TRUE,
                                  optimize_cutoffs = TRUE,
                                  p = 0.015,
                                  lfc = 1.5,
                                  gene_expression_threshold = 0.5) {
  if (is.null(single_cell_obj)) {
    stop("Parameter 'single_cell_obj' is missing or null, but is required.")
  }

  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is missing or null, but is required.")
  }

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

  if (!assay_sp %in% names(SummarizedExperiment::assays(spatial_obj))) {
    message(
      "requested assay ",
      assay_sp,
      " not available in spatial object. Using first available assay."
    )
    assay_sp <- names(SummarizedExperiment::assays(spatial_obj))[1]
  }

  ad <- spe_to_ad(single_cell_obj, assay = assay_sc)
  ad$obs[[cell_type_col]] <- as.character(SingleCellExperiment::colData(single_cell_obj)[[cell_type_col]])
  ad$obs_names <- colnames(single_cell_obj)
  ad$var_names <- rownames(single_cell_obj)

  bulk_matrix <- Matrix::t(SummarizedExperiment::assay(spatial_obj, assay_sp))
  bulk_df <- as.data.frame(as.matrix(bulk_matrix))
  rownames(bulk_df) <- colnames(spatial_obj)
  colnames(bulk_df) <- rownames(spatial_obj)

  reticulate::source_python(system.file("python", "rectangle_spacedeconv.py", package = "spacedeconv"))

  deconv <- py_deconvolute_rectangle(
    sc_obj = ad,
    bulks = bulk_df,
    cell_type_col = cell_type_col,
    layer = layer,
    raw = raw,
    correct_mrna_bias = correct_mrna_bias,
    optimize_cutoffs = optimize_cutoffs,
    p = p,
    lfc = lfc,
    n_cpus = n_cpus,
    gene_expression_threshold = gene_expression_threshold
  )

  deconv <- as.matrix(deconv)
  if (nrow(deconv) == ncol(spatial_obj)) {
    rownames(deconv) <- colnames(spatial_obj)
  }
  deconv <- attachToken(deconv, result_name)

  return(deconv)
}
