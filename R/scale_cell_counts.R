#' Scale Relative Cell Fractions to Absolute Counts
#'
#' Multiplies a relative fraction column by a cell-count column and stores the
#' result in `colData`.
#'
#' @param spe A `SpatialExperiment`.
#' @param value Column name in `colData(spe)` with relative fractions.
#' @param cell_counts Column name in `colData(spe)` with absolute cell counts.
#' @param resName Optional name for the output column. Defaults to
#' `paste0(value, "_absolute")`.
#'
#' @return `SpatialExperiment` with the scaled counts added to `colData`.
#'
#' @export
scale_cell_counts <- function(spe, value, cell_counts, resName = NULL) {
  # Validate inputs
  if (!value %in% colnames(colData(spe))) {
    stop("Specified value column does not exist in colData of the SPE object.")
  }

  if (is.null(resName)) {
    resName <- paste0(value, "_absolute")
  }

  if (cell_counts %in% colnames(colData(spe))) {
    cell_densities <- colData(spe)[, cell_counts]
  } else {
    stop("Cell count data not available in the object")
  }

  # Calculate absolute cell densities
  absoluteValues <- colData(spe)[, value] * cell_densities

  # Add the absolute cell densities as a new column
  colData(spe)[[resName]] <- absoluteValues

  return(spe)
}
