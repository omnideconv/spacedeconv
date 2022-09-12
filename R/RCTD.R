#' RCTD builds the model internally, please just use the deconvolute method
#'
#' @return NULL
#' @export
build_model_rctd <- function() {
  message("RCTD builds it's model internally, please just use the deconvolute method")
  return(NULL)
}


#' RCTD Deconvolution
#' @param sce SingleCellExperiment
#' @param cell_type_col Column containting cell type annotation
#' @param spe SpatialExperiment
#' @param assay_sc single cell assay to use
#' @param assay_sp spatial assay to use
#' @param n_umi_sc (optional) named list of umi counts for each cell
#' @param n_umi_sp (optional) named list of umi counts for each spot
#' @param n_cores Number of CPU cores to use for the calculation, NULL = use all cores
deconvolute_rctd <- function(sce, cell_type_col, spe, assay_sc = "counts", assay_sp = "counts", n_umi_sc = NULL, n_umi_sp = NULL, n_cores = NULL) {
  if (is.null(sce)) {
    stop("Parameter 'sce' is missing or null, but is required.")
  }

  if (is.null(cell_type_col)) {
    stop("Parameter 'cell_type_col' is missing or null, but is required.")
  }

  if (!is.null(n_cores) && n_cores < 1) {
    stop("Parameter 'n_cores' needs to be a positive integer")
  }

  if (!checkCol(sce, cell_type_col)) {
    stop(paste0("Column \"", cell_type_col, "\" can't be found in single cell object"))
  }

  message("Preparing Data for RCTD")

  # check if requested assay exists
  if (!assay_sc %in% names(SummarizedExperiment::assays(sce))) {
    message(
      "requested assay ",
      assay_sc,
      " not available in expression object. Using first available assay."
    )
    assay_sc <- names(SummarizedExperiment::assays(sce))[1] # change to first available assay request not available
  }

  # check if requested assay exists
  if (!assay_sp %in% names(SummarizedExperiment::assays(spe))) {
    message(
      "requested assay ",
      assay_sp,
      " not available in expression object. Using first available assay."
    )
    assay_sp <- names(SummarizedExperiment::assays(spe))[1] # change to first available assay request not available
  }

  # sc counts
  sc_counts <- SummarizedExperiment::assay(sce, assay_sc)
  sc_counts <- methods::as(sc_counts, "dgCMatrix") # ensure datatype
  # in this specific case we must make gene names unique! There are duplicates


  # named cell type factor
  cell_types <- SingleCellExperiment::colData(sce)[, cell_type_col] # get cell type factor
  names(cell_types) <- rownames(colData(sce)) # name it
  cell_types <- as.factor(cell_types) # ensure factor

  # spatial_counts
  spatial_counts <- SummarizedExperiment::assay(spe, assay_sp)
  spatial_counts <- methods::as(spatial_counts, "dgCMatrix")
  # in this specific case we must make gene names unique! There are duplicates
  rownames(spatial_counts) <- make.names(rownames(spatial_counts), unique = T)

  # spatial coords
  spatial_coords <- SpatialExperiment::spatialCoords(spe)
  spatial_coords <- as.data.frame(spatial_coords) # ensure dataframe for RCTD

  # n_umi_sc Optional!!!
  # n_umi_sp, optional!!!


  # n_cores, if null = max available
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores()
  }


  # create reference object
  reference <- spacexr::Reference(
    counts = sc_counts,
    cell_types = cell_types,
    nUMI = n_umi_sc
  )

  # create spatial dataset
  puck <- spacexr::SpatialRNA(
    coords = spatial_coords,
    counts = spatial_counts,
    nUMI = n_umi_sp
  )

  # create rctd dataset
  rctd_object <- spacexr::create.RCTD(
    spatialRNA = puck,
    reference = reference,
    max_cores = n_cores,
    CELL_MIN_INSTANCE = 1
    # UMI_min = 0
  )

  message("Starting RCTD Deconvolution")

  # perform deconvolution
  rctd_object <- spacexr::run.RCTD(rctd_object)

  # manage results
  results <- rctd_object@results
  normalized_results <- spacexr::normalize_weights(results$weights)
  normalized_results <- as.matrix(normalized_results)

  return(normalized_results)
}


# #' RCTD Deconvolution
# #' @param sc_counts (dgC) matrix of single cell counts
# #' @param cell_types Named factor of cell types for the sc_counts matrix
# #' @param n_umi_sc (optional) named list of umi counts for each cell
# #' @param spatial_counts (dgC) matrix of spatial counts
# #' @param spatial_coords data frame containing spatial coordinates
# #' @param n_umi_sp (optional) named list of umi counts for each spot
# #' @param n_cores Number of CPU cores to use for the calculation
# #' @export
# # deconvolute_rctd <- function(sc_counts, cell_types, n_umi_sc = NULL, spatial_counts, spatial_coords, n_umi_sp = NULL, n_cores = 20) {
#   if (is.null(sc_counts)) {
#     stop("Parameter 'sc_counts' is missing or null, but is required.")
#   }
#
#   if (is.null(cell_types)) {
#     stop("Parameter 'cell_types' is missing or null, but is required.")
#   }
#
#   if (is.null(spatial_counts)) {
#     stop("Parameter 'spatial_counts' is missing or null, but is required.")
#   }
#
#   if (is.null(spatial_coords)) {
#     stop("Parameter 'spatial_coords' is missing or null, but is required.")
#   }
#   if (n_cores < 1) {
#     stop("Parameter 'n_cores' needs to be a positive integer")
#   }
#
#   # create reference object
#   reference <- spacexr::Reference(
#     counts = sc_counts,
#     cell_types = cell_types,
#     nUMI = n_umi_sc
#   )
#
#   # create spatial dataset
#   puck <- spacexr::SpatialRNA(
#     coords = spatial_coords,
#     counts = spatial_counts,
#     nUMI = n_umi_sp
#   )
#
#   # create rctd dataset
#   rctd_object <- spacexr::create.RCTD(
#     spatialRNA = puck,
#     reference = reference,
#     max_cores = parallel::detectCores() - 1
#   )
#
#   # perform deconvolution
#   rctd_object <- spacexr::run.RCTD(rctd_object)
#
#   # manage results
#   results <- rctd_object@results
#   normalized_results <- spacexr::normalize_weights(results$weights)
#   normalized_results <- as.matrix(normalized_results)
#
#   return(normalized_results)
# }
