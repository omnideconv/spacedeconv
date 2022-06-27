#' RCTD Deconvolution
#' @param sc_counts (dgC) matrix of single cell counts
#' @param cell_types Named factor of cell types for the sc_counts matrix
#' @param n_umi_sc (optional) named list of umi counts for each cell
#' @param spatial_counts (dgC) matrix of spatial counts
#' @param spatial_coords data frame containing spatial coordinates
#' @param n_umi_sp (optional) named list of umi counts for each spot
#' @param n_cores Number of CPU cores to use for the calculation
#' @export

deconvolute_rctd <- function(sc_counts, cell_types, n_umi_sc = NULL, spatial_counts, spatial_coords, n_umi_sp = NULL, n_cores = 20) {
  if (is.null(sc_counts)) {
    stop("Parameter 'sc_counts' is missing or null, but is required.")
  }

  if (is.null(cell_types)) {
    stop("Parameter 'cell_types' is missing or null, but is required.")
  }

  if (is.null(spatial_counts)) {
    stop("Parameter 'spatial_counts' is missing or null, but is required.")
  }

  if (is.null(spatial_coords)) {
    stop("Parameter 'spatial_coords' is missing or null, but is required.")
  }
  if (n_cores < 1) {
    stop("Parameter 'n_cores' needs to be a positive integer")
  }

  # create reference
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
    max_cores = parallel::detectCores() - 1
  )

  # perform deconvolution
  rctd_object <- spacexr::run.RCTD(rctd_object)

  # manage results
  results <- rctd_object@results
  normalized_results <- spacexr::normalize_weights(results$weights)
  normalized_results <- as.matrix(normalized_results)

  return(normalized_results)
}
