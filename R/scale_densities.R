#' Scale relative cell fractions to absolute by providing cell densities
#'
#' This function calculates absolute cell densities for a SpatialExperiment object
#' by scaling a specified column of relative cell fractions with provided absolute cell densities.
#' The result is added as a new column in the `colData` of the SPE object.
#'
#' @param spe An SpatialExperiment object containing spatial transcriptomics data.
#' @param value A string specifying the name of the column in `colData(spe)` that contains
#'              relative cell fractions/deconvolution results to be scaled.
#' @param cell_densities A named vector of absolute cell densities, where names correspond
#'                       to the rownames in `colData(spe)`.
#' @param resName Optional; a string specifying the name of the new column to be added
#'                to `colData(spe)` that will contain the absolute cell densities.
#'                If NULL, defaults to appending "_absolute" to the `value` parameter.
#'
#' @return SPE object with added absolute celltype counts for each spot
#'
#' @export
scale_cell_densities <- function(spe, value, cell_densities, resName = NULL) {
  # Validate inputs
  if (!value %in% colnames(colData(spe))) {
    stop("Specified value column does not exist in colData of the SPE object.")
  }

  if (is.null(resName)) {
    resName <- paste0(value, "_absolute")
  }

  if (any(is.na(match(rownames(colData(spe)), names(cell_densities))))) {
    stop("Not all rownames in colData have a matching entry in cell_densities.")
  }

  # Calculate absolute cell densities
  absoluteValues <- colData(spe)[, value] * cell_densities[rownames(colData(spe))]

  # Add the absolute cell densities as a new column
  colData(spe)[[resName]] <- absoluteValues

  return(spe)
}
