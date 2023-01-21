#' Aggregate Deconvolution results
#'
#' Aggregate deconvolution results by merging result columns in the spatial object
#'
#' @param spatial_obj SpatialExperiment containing deconvolution results
#' @param cell_type_1 cell type to aggregate, including the method "method_celltype"
#' @param cell_type_2 cell type to aggregate, including the method "method_celltype"
#' @param name new name for aggregation
#' @param remove logical, remove provided cell types and just keep the aggregation
#'
#' @returns SpatialObject containing aggregation of provided cell types
aggregate <- function(spatial_obj = NULL, cell_type_1 = NULL, cell_type_2 = NULL,
                      name = NULL, remove = FALSE) {
  cli::cli_rule(left = "spacedeconv")

  cli::cli_progress_step("testing parameter", msg_done = "parameter OK")

  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is null or missing, but is required")
  }

  if (is.null(cell_type_1) || is.null(cell_type_2)) {
    stop("One celltype is NULL, please provide two celltypes")
  }

  # check if cell types are available in object
  if (!cell_type_1 %in% names(colData(spatial_obj))) {
    stop(paste("celltype", cell_type_1, "is not available in spatial object"))
  }

  if (!cell_type_2 %in% names(colData(spatial_obj))) {
    stop(paste("celltype", cell_type_2, "is not available in spatial object"))
  }

  # check if cell types are different
  if (cell_type_1 == cell_type_2) {
    stop("provided cell types need to be distinct!")
  }

  # check and calculate new name
  if (is.null(name)) {
    name <- paste0(cell_type_1, cell_type_2)
  }

  cli::cli_progress_step("Aggregating cell types", msg_done = "Aggregated cell types")

  # aggregate and set colnames
  aggregation <- colData(spatial_obj)[, cell_type_1] + colData(spatial_obj)[, cell_type_2]

  # remove old cell types if requested
  if (remove) {
    colData(spatial_obj) <- colData(spatial_obj)[, !names(colData(spatial_obj)) %in% c(cell_type_1, cell_type_2)]
  }

  # merge data
  tmp <- cbind(colData(spatial_obj), aggregation)
  colnames(tmp) <- c(names(colData(spatial_obj)), name)
  colData(spatial_obj) <- tmp

  cli::cli_process_done()

  return(spatial_obj)
}
