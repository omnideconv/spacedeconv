#' Build a SPOTlight Model
#'
#' Trains the SPOTlight model using single-cell and spatial data. Marker
#' genes must be supplied explicitly; automatic marker discovery is not supported.
#'
#' @param single_cell_obj `SingleCellExperiment`.
#' @param cell_type_col Column with cell type labels.
#' @param spatial_obj `SpatialExperiment`.
#' @param assay_sc Single-cell assay to use.
#' @param assay_sp Spatial assay to use.
#' @param markers Required data.frame or S4Vectors DataFrame with columns `gene`,
#' `cluster` and numeric `mean.AUC` weights in (0, 1]. Cluster labels must match
#' `cell_type_col`. Each cell type needs a marker present in both expression
#' objects. Gene/cluster pairs must be unique. `NULL` raises an error.
#' @param ... Additional training parameters passed to `SPOTlight::trainNMF()`.
build_model_spotlight <- function(single_cell_obj, cell_type_col = "cell_ontology_class", spatial_obj, assay_sc = "counts", assay_sp = "counts", markers = NULL, ...) {
  if (is.null(single_cell_obj)) {
    stop("Parameter 'single_cell_obj' is null or missing, but is required")
  }

  if (!checkCol(single_cell_obj, cell_type_col)) {
    stop(paste0("Column \"", cell_type_col, "\" can't be found in single cell object"))
  }

  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial'obj' is null or missing, but is required")
  }

  # check if requested assay exists
  if (!assay_sc %in% names(SummarizedExperiment::assays(single_cell_obj))) {
    message(
      "requested assay ",
      assay_sc,
      " not available in expression object. Using first available assay."
    )
    assay_sc <- names(SummarizedExperiment::assays(single_cell_obj))[1] # change to first available assay request not available
  }

  # check if requested assay exists
  if (!assay_sp %in% names(SummarizedExperiment::assays(spatial_obj))) {
    message(
      "requested assay ",
      assay_sp,
      " not available in expression object. Using first available assay."
    )
    assay_sp <- names(SummarizedExperiment::assays(spatial_obj))[1] # change to first available assay request not available
  }

  groups <- as.character(SummarizedExperiment::colData(single_cell_obj)[[cell_type_col]])
  if (is.null(markers)) {
    stop("SPOTlight requires 'markers': supply a data frame with gene, cluster and mean.AUC columns. Automatic marker discovery is not supported.", call. = FALSE)
  }
  if (!(is.data.frame(markers) || methods::is(markers, "DataFrame")) ||
      !all(c("gene", "cluster", "mean.AUC") %in% names(markers)) || nrow(markers) == 0L) {
    stop("SPOTlight 'markers' must be a non-empty data frame with gene, cluster and mean.AUC columns.", call. = FALSE)
  }
  mgs <- as.data.frame(markers)
  for (column in c("gene", "cluster")) {
    if (!(is.character(mgs[[column]]) || is.factor(mgs[[column]])) ||
        anyNA(mgs[[column]]) || any(!nzchar(trimws(as.character(mgs[[column]]))))) {
      stop("SPOTlight marker gene and cluster labels must be non-empty strings without NA.", call. = FALSE)
    }
    mgs[[column]] <- as.character(mgs[[column]])
  }
  if (!is.numeric(mgs$mean.AUC) || any(!is.finite(mgs$mean.AUC)) ||
      any(mgs$mean.AUC <= 0 | mgs$mean.AUC > 1)) {
    stop("SPOTlight marker mean.AUC weights must be finite numeric values in (0, 1].", call. = FALSE)
  }
  if (anyDuplicated(mgs[c("gene", "cluster")])) {
    stop("SPOTlight markers contain duplicate gene/cluster pairs.", call. = FALSE)
  }
  if (anyNA(groups) || any(!nzchar(trimws(groups))) ||
      !setequal(unique(mgs$cluster), unique(groups))) {
    stop("SPOTlight marker clusters must match all cell types in 'cell_type_col', without missing labels.", call. = FALSE)
  }
  shared_genes <- intersect(rownames(single_cell_obj), rownames(spatial_obj))
  covered <- mgs$cluster[mgs$gene %in% shared_genes]
  if (!all(unique(groups) %in% covered)) {
    stop("SPOTlight requires at least one marker per cell type present in both expression objects.", call. = FALSE)
  }

  model <- SPOTlight::trainNMF(
    x = single_cell_obj,
    y = spatial_obj,
    groups = groups,
    mgs = mgs,
    weight_id = "mean.AUC",
    slot_sc = assay_sc, # not sure about this one!
    slot_sp = assay_sp,
    ...
  )

  return(model)
}


#' Deconvolute with SPOTlight
#'
#' Runs SPOTlight deconvolution and returns a deconvolution matrix with column
#' names prefixed by `result_name`.
#'
#' @param spatial_obj `SpatialExperiment`.
#' @param model SPOTlight model from `build_model_spotlight()`.
#' @param assay_sp Spatial assay to use.
#' @param result_name Prefix used to label result columns (default: "spotlight").
deconvolute_spotlight <- function(spatial_obj, model = NULL, assay_sp = "counts", result_name = "spotlight") {
  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' missing or null, but is required")
  }

  if (is.null(model)) {
    stop("Model is missing or null, but is required")
  }

  # check for model integrity: names(signature) must be "mod" "topic"

  # extract model information
  mod <- model$mod
  ref <- model$topic

  # SPOTlight 1.6.3 passes its cell-type-by-topic reference directly to NNLS,
  # although NNLS needs topics as rows and cell types as columns. Its output
  # labels also come from reference rows rather than coefficient columns.
  cell_types <- rownames(ref)
  topics <- rownames(NMF::coef(mod))
  if (is.null(cell_types) || is.null(topics) ||
      anyDuplicated(cell_types) || !setequal(colnames(ref), topics)) {
    stop("SPOTlight reference must have unique cell-type names and match the model topics.",
         call. = FALSE)
  }
  ref <- t(ref[, topics, drop = FALSE])

  # deconvolute
  deconv <- SPOTlight::runDeconvolution(
    x = spatial_obj,
    mod = mod,
    ref = ref,
    slot = assay_sp
  )

  deconvolution <- deconv$mat
  colnames(deconvolution) <- cell_types

  # attach method token
  deconvolution <- attachToken(deconvolution, result_name)

  return(deconvolution)
}
