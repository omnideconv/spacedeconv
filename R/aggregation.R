#' Aggregate Deconvolution Results
#'
#' Aggregates deconvolution results by merging result columns in the spatial object.
#' This function now supports aggregating multiple cell types by passing a list of
#' cell types. The older method of specifying two separate cell types is deprecated
#' but still supported for backward compatibility.
#'
#' @param spatial_obj SpatialExperiment containing deconvolution results.
#' @param cell_types A list of cell types to aggregate. This is the preferred method
#' for specifying cell types to aggregate. If this parameter is used, 'cell_type_1'
#' and 'cell_type_2' are ignored.
#' @param cell_type_1 (Deprecated) The first cell type to aggregate, included for
#' backward compatibility. Use 'cell_types' for new code.
#' @param cell_type_2 (Deprecated) The second cell type to aggregate, included for
#' backward compatibility. Use 'cell_types' for new code.
#' @param name Optional new name for the aggregated result. If not provided, a name
#' is generated automatically by concatenating the cell type names.
#' @param remove Logical, specifies whether to remove the provided cell types from
#' the spatial object after aggregation.
#'
#' @return A SpatialExperiment object containing the aggregation of the provided
#' cell types.
#'
#' @details The function will aggregate the specified cell types by summing their
#' corresponding columns in the spatial object. If the deprecated parameters
#' 'cell_type_1' and 'cell_type_2' are used, a warning is issued, and these types
#' are added to the 'cell_types' list for aggregation.
#'
#' @examples
#' \dontrun{
#' # For new code, use:
#' aggregate_results(spatial_data_3, cell_types = c("cell_type_A", "cell_type_B"))
#'
#' # For backward compatibility:
#' aggregate_results(spatial_data_3, cell_type_1 = "cell_type_A", cell_type_2 = "cell_type_B")
#' }
#'
#' @export
aggregate_results <- function(spatial_obj = NULL, cell_types = NULL, cell_type_1 = NULL, cell_type_2 = NULL,
                              name = NULL, remove = FALSE) {
  cli::cli_rule(left = "spacedeconv")

  cli::cli_progress_step("testing parameter", msg_done = "parameter OK")

  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is null or missing, but is required")
  }

  # Handle deprecated parameters
  if (!is.null(cell_type_1) || !is.null(cell_type_2)) {
    warning("Parameters 'cell_type_1' and 'cell_type_2' are deprecated. Please use 'cell_types' instead.")
    cell_types <- unique(c(cell_types, cell_type_1, cell_type_2))
    cell_types <- cell_types[!is.null(cell_types)] # Remove NULLs
  }

  if (is.null(cell_types) || length(cell_types) < 2) {
    stop("Please provide a list of at least two cell types")
  }

  # check if cell types are available in object
  unavailable_types <- cell_types[!cell_types %in% names(colData(spatial_obj))]
  if (length(unavailable_types) > 0) {
    stop(paste("cell types", paste(unavailable_types, collapse = ", "), "are not available in spatial object"))
  }

  # check for duplicate cell types
  if (length(unique(cell_types)) != length(cell_types)) {
    stop("Duplicate cell types found in the list")
  }

  # check and calculate new name
  if (is.null(name)) {
    name <- paste(cell_types, collapse = "_")
  }

  cli::cli_progress_step("Aggregating cell types", msg_done = "Aggregated cell types")

  # aggregate cell types
  aggregation <- rowSums(as.matrix(colData(spatial_obj)[, cell_types]))

  # remove old cell types if requested
  if (remove) {
    colData(spatial_obj) <- colData(spatial_obj)[, !names(colData(spatial_obj)) %in% cell_types]
  }

  # set new name
  aggregation <- as.data.frame(aggregation)
  colnames(aggregation) <- name

  # merge data
  tmp <- cbind(colData(spatial_obj), aggregation)
  colData(spatial_obj) <- tmp

  cli::cli_process_done()

  return(spatial_obj)
}
