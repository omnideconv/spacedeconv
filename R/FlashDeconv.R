#' FlashDeconv: Fast Linear Algebra for Scalable Hybrid Deconvolution
#'
#' FlashDeconv is a high-performance spatial transcriptomics deconvolution method
#' designed for atlas-scale data. It uses structure-preserving randomized sketching
#' and spatial graph regularization for efficient deconvolution.
#'
#' Key features:
#' - Ultra-fast: Processes 1 million spots in ~3 minutes
#' - Linear O(N) scaling for time and memory
#' - No GPU required
#' - Rare cell detection via leverage-score sampling
#'
#' @references
#' Yang, C., Chen, J. & Zhang, X. FlashDeconv enables atlas-scale, multi-resolution
#' spatial deconvolution via structure-preserving sketching. bioRxiv (2025).
#' https://doi.org/10.64898/2025.12.22.696108
#'
#' @seealso
#' GitHub: https://github.com/cafferychen777/flashdeconv
#'
#' @name flashdeconv
NULL


#' No signature building required for FlashDeconv
#'
#' FlashDeconv performs deconvolution directly without requiring a pre-built
#' signature matrix. Just call the deconvolute method directly.
#'
#' @returns NULL
#' @export
build_model_flashdeconv <- function() {
  message("FlashDeconv does not require pre-building a signature. Just call the deconvolute method directly.")
  return(NULL)
}


#' Deconvolute using FlashDeconv
#'
#' Perform cell type deconvolution on spatial transcriptomics data using FlashDeconv.
#'
#' @param single_cell_obj SingleCellExperiment containing reference single-cell expression data
#' @param spatial_obj SpatialExperiment containing spatial transcriptomics data
#' @param cell_type_col Column name in single_cell_obj containing cell type annotations
#' @param assay_sc Assay name to use from single_cell_obj, default = "counts"
#' @param assay_sp Assay name to use from spatial_obj, default = "counts"
#' @param sketch_dim Dimension of the sketched space. Higher values preserve more information
#'   but increase computation. Default: 512
#' @param lambda_spatial Spatial regularization strength. Higher values encourage smoother
#'   spatial patterns. Default: 5000. Typical range: 1000-10000 for Visium.
#' @param n_hvg Number of highly variable genes to select. Default: 2000
#' @param n_markers_per_type Number of marker genes per cell type. Default: 50
#' @param result_name Token to identify deconvolution results in object, default = "flashdeconv"
#'
#' @returns A data.frame with cell type proportions (spots x cell types)
#'
#' @export
deconvolute_flashdeconv <- function(single_cell_obj,
                                    spatial_obj,
                                    cell_type_col = "cell_ontology_class",
                                    assay_sc = "counts",
                                    assay_sp = "counts",
                                    sketch_dim = 512L,
                                    lambda_spatial = 5000.0,
                                    n_hvg = 2000L,
                                    n_markers_per_type = 50L,
                                    result_name = "flashdeconv") {
  # Parameter validation
  if (is.null(single_cell_obj)) {
    stop("Parameter 'single_cell_obj' is missing or null, but is required!")
  }

  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is missing or null, but is required!")
  }

  if (!cell_type_col %in% names(colData(single_cell_obj))) {
    stop(paste0("Column '", cell_type_col, "' not found in single_cell_obj"))
  }

  # Check if requested assays exist
  if (!assay_sc %in% names(SummarizedExperiment::assays(single_cell_obj))) {
    message(
      "Requested assay '", assay_sc,
      "' not available in single_cell_obj. Using first available assay."
    )
    assay_sc <- names(SummarizedExperiment::assays(single_cell_obj))[1]
  }

  if (!assay_sp %in% names(SummarizedExperiment::assays(spatial_obj))) {
    message(
      "Requested assay '", assay_sp,
      "' not available in spatial_obj. Using first available assay."
    )
    assay_sp <- names(SummarizedExperiment::assays(spatial_obj))[1]
  }

  # Convert spatial object to AnnData
  sp_ad <- spe_to_ad_flashdeconv(spatial_obj, assay = assay_sp)

  # Convert single cell object to AnnData
  ref_ad <- sce_to_ad_flashdeconv(single_cell_obj, assay = assay_sc, cell_type_col = cell_type_col)

  # Source the Python script
  reticulate::source_python(
    system.file("python", "flashdeconv_spacedeconv.py", package = "spacedeconv")
  )

  # Call the Python deconvolution function
  deconv <- py_deconvolute_flashdeconv(
    sp_obj = sp_ad,
    ref_obj = ref_ad,
    cell_type_col = cell_type_col,
    sketch_dim = as.integer(sketch_dim),
    lambda_spatial = as.numeric(lambda_spatial),
    n_hvg = as.integer(n_hvg),
    n_markers_per_type = as.integer(n_markers_per_type),
    layer_st = NULL,
    layer_ref = NULL,
    spatial_key = "spatial"
  )

  # Make column names valid R names
  colnames(deconv) <- make.names(colnames(deconv))

  # Attach result token
  deconv <- attachToken(deconv, result_name)

  return(deconv)
}


#' Convert SpatialExperiment to AnnData for FlashDeconv
#'
#' Internal function to convert SpatialExperiment to AnnData format
#' with spatial coordinates properly set.
#'
#' @param spe SpatialExperiment object
#' @param assay Assay to use
#'
#' @returns AnnData object with spatial coordinates
#' @keywords internal
spe_to_ad_flashdeconv <- function(spe, assay = "counts") {
  if (is.null(spe)) {
    stop("Parameter 'spe' is missing or null, but is required")
  }

  if (!assay %in% names(SummarizedExperiment::assays(spe))) {
    stop("Requested assay not available in object")
  }

  # Get expression matrix (genes x spots -> transpose to spots x genes for AnnData)
  expr_matrix <- Matrix::t(SummarizedExperiment::assay(spe, assay))

  # Get spatial coordinates
  spatial_coords <- SpatialExperiment::spatialCoords(spe)

  # Create AnnData object
  ad <- anndata::AnnData(
    X = expr_matrix,
    var = as.data.frame(SingleCellExperiment::rowData(spe)),
    obs = as.data.frame(SingleCellExperiment::colData(spe))
  )

  # Add spatial coordinates to obsm
  ad$obsm[["spatial"]] <- spatial_coords

  return(ad)
}


#' Convert SingleCellExperiment to AnnData for FlashDeconv
#'
#' Internal function to convert SingleCellExperiment to AnnData format.
#'
#' @param sce SingleCellExperiment object
#' @param assay Assay to use
#' @param cell_type_col Column containing cell type annotations
#'
#' @returns AnnData object
#' @keywords internal
sce_to_ad_flashdeconv <- function(sce, assay = "counts", cell_type_col = "cell_type") {
  if (is.null(sce)) {
    stop("Parameter 'sce' is missing or null, but is required")
  }

  if (!assay %in% names(SummarizedExperiment::assays(sce))) {
    stop("Requested assay not available in object")
  }

  # Get expression matrix (genes x cells -> transpose to cells x genes for AnnData)
  expr_matrix <- Matrix::t(SummarizedExperiment::assay(sce, assay))

  # Get metadata
  obs_data <- as.data.frame(SingleCellExperiment::colData(sce))
  var_data <- as.data.frame(SingleCellExperiment::rowData(sce))

  # Ensure cell type column is available
  if (!cell_type_col %in% names(obs_data)) {
    stop(paste0("Cell type column '", cell_type_col, "' not found in object"))
  }

  # Create AnnData object
  ad <- anndata::AnnData(
    X = expr_matrix,
    var = var_data,
    obs = obs_data
  )

  return(ad)
}
