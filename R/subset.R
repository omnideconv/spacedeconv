#' Subset a SingleCellExperiment by Cell Type
#'
#' Randomly subsamples cells from a `SingleCellExperiment` either evenly across
#' cell types or proportionally to their original abundance. If a cell type has
#' fewer cells than requested, you can keep all of them (`"asis"`) or drop that
#' cell type (`"remove"`).
#'
#' @param sce A `SingleCellExperiment` with cell type annotations in `colData`.
#' @param cell_type_col Column name in `colData` containing cell type labels.
#' @param scenario Subsampling strategy: `"even"` selects equal counts per type,
#' `"mirror"` preserves original proportions.
#' @param ncells Target total number of cells to keep (default: 1000).
#' @param notEnough Behavior when a cell type has too few cells: `"asis"` keeps all,
#' `"remove"` drops that cell type.
#' @param seed Random seed for reproducible sampling.
#'
#' @return A subsetted `SingleCellExperiment`.
#'
#' @export
subsetSCE <- function(sce, cell_type_col = "celltype_major", scenario = "even", ncells = 1000, notEnough = "asis", seed = 12345) {
  cli::cli_rule(left = "spacedeconv")

  cli::cli_progress_step("testing parameter", msg_done = "parameter OK")

  if (is.null(sce)) {
    stop("SingleCellExperiment missing but required")
  }

  # check that row and colnames are set
  if (checkRowColumn(sce)) {
    stop("Rownames or colnames are not set for SingleCellExperiment but are required")
  }

  # check that cell_type_col exists in object
  if (!checkCol(sce, cell_type_col)) {
    stop(paste0("Column \"", cell_type_col, "\" can't be found in single cell object"))
  }

  if (!is.numeric(seed)) {
    stop("seed has to be numeric!")
  }

  set.seed(seed)
  cli::cli_alert_info(paste0("Set seed to ", seed))

  cli::cli_progress_step(msg = paste0("extracting up to ", ncells, " cells"))

  # full selection vector
  x <- rep(FALSE, ncol(sce))

  # scenarios, build selection vector
  # general idea:
  # get the cells for each celltype, select a suitable number
  # and set these positions in x to TRUE, then
  # subset with x
  if (scenario == "even") {
    # calculate number of cells per celltype
    cells_per_type <- floor(ncells / length(unique(sce[[cell_type_col]])))

    for (celltype in unique(sce[[cell_type_col]])) {
      locations <- which(sce[[cell_type_col]] == celltype)

      # check if enough cells available
      # if not use all
      if (length(locations) >= cells_per_type) {
        # select "cells_per_type" random cells
        x[sample(locations, size = cells_per_type, replace = FALSE)] <- TRUE
      } else {
        if (notEnough == "asis") {
          # just use all
          cli::cli_alert_info(paste0("Not enough cells for ", celltype, ". Using all avialable cells"))

          # use all cells
          x[locations] <- TRUE
        } else if (notEnough == "remove") {
          cli::cli_alert_warning(paste0("Not enough cells for ", celltype, ". Removing Celltype"))
          x[locations] <- FALSE
        }
      }
    }
  } else if (scenario == "mirror") {
    # Calculate total counts per cell type
    cell_counts <- table(sce[[cell_type_col]])

    # Calculate proportions
    cell_proportions <- cell_counts / sum(cell_counts)

    # Calculate the number of cells for each cell type based on proportions
    cells_per_type <- round(cell_proportions * ncells)

    for (celltype in names(cells_per_type)) {
      locations <- which(sce[[cell_type_col]] == celltype)
      num_cells_to_select <- cells_per_type[celltype]

      if (length(locations) >= num_cells_to_select) {
        x[sample(locations, size = num_cells_to_select, replace = FALSE)] <- TRUE
      } else {
        if (notEnough == "asis") {
          cli::cli_alert_info(paste0("Not enough cells for ", celltype, ". Using all available cells"))
          x[locations] <- TRUE
        } else if (notEnough == "remove") {
          cli::cli_alert_warning(paste0("Not enough cells for ", celltype, ". Removing Celltype"))
          x[locations] <- FALSE
        }
      }
    }
  }


  # actually subset and return the sce
  sce <- sce[, x]
  cli::cli_progress_step(paste("extracted", ncol(sce), "cells"))

  cli::cli_progress_done()

  return(sce)
}


#' Subset a SpatialExperiment by Spatial Coordinates
#'
#' Filters a `SpatialExperiment` to spots that fall within the provided x/y ranges
#' of the spatial coordinates. If a range is `NULL`, the full range is used.
#'
#' @param spe A `SpatialExperiment` object.
#' @param colRange Numeric vector of length 2 giving min/max x-coordinates
#' (`pxl_col_in_fullres`).
#' @param rowRange Numeric vector of length 2 giving min/max y-coordinates
#' (`pxl_row_in_fullres`).
#'
#' @return A `SpatialExperiment` subsetted to the specified coordinate window.
#'
#' @export
subsetSPE <- function(spe, colRange = NULL, rowRange = NULL) {
  # Extract spatial coordinates
  coords <- spatialCoords(spe)

  # Use full range if colRange or rowRange is NULL, sort if in wrong order
  if (is.null(colRange)) {
    colRange <- range(coords[, 1], na.rm = TRUE)
  } else {
    colRange <- sort(colRange)
  }

  if (is.null(rowRange)) {
    rowRange <- range(coords[, 2], na.rm = TRUE)
  } else {
    rowRange <- sort(rowRange)
  }

  # Subset the coordinates based on the specified range
  subsetCoords <- coords[coords[, 1] >= colRange[1] & coords[, 1] <= colRange[2] &
    coords[, 2] >= rowRange[1] & coords[, 2] <= rowRange[2], ]

  # Find the indices of the subset
  indices <- match(rownames(subsetCoords), rownames(coords))

  # Subset the SpatialExperiment object
  fovSubset <- spe[, indices]

  return(fovSubset)
}
