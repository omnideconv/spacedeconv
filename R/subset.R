#' Subset a SingleCellExperiment by Cell Type
#'
#' This function facilitates the subsetting of a `SingleCellExperiment` object to a specified number of cells,
#' with an option to maintain an even or proportional distribution across different cell types. It offers flexibility
#' in handling scenarios where insufficient cells of a specific type are available and ensures reproducibility through
#' the option to set a seed for random sampling.
#'
#' @param sce A `SingleCellExperiment` object containing cell type annotations. The function applies subsetting
#' operations on this object to achieve a targeted composition of cell types.
#' @param cell_type_col A character string indicating the column name within `colData` of the `SingleCellExperiment`
#' object that contains cell type information. This column is used to identify and categorize cells by type for subsetting.
#' @param scenario Specifies the strategy for cell selection across types: `"even"` for equal numbers from each cell
#' type, or `"mirror"` to reflect the proportional representation of cell types in the subset as in the original dataset.
#' @param ncells The target total number of cells for the subsetted object. Defaults to 1000. This parameter dictates
#' the overall size of the subsetted `SingleCellExperiment`.
#' @param notEnough Defines the behavior when there are not enough cells of a certain type to meet the subsetting criteria:
#' `"remove"` to exclude that cell type entirely, or `"asis"` to include all available cells of that type despite falling short
#' of the target number.
#' @param seed An optional numeric value to set the random seed for reproducibility of the subsetting process. This ensures
#' that the subset selection can be replicated in future analyses.
#'
#' @return Returns a subsetted `SingleCellExperiment` object. The resulting object reflects the targeted distribution and
#' number of cells as specified by the input parameters. It enables focused analysis on a more manageable dataset while
#' preserving the desired representation of cell types.
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


#' Isolate Regions within a SpatialExperiment Object
#'
#' Tailor a `SpatialExperiment` object to focus on specific areas by defining ranges in the spatial coordinates.
#' This function narrows down the dataset to include only those spots that fall within user-defined x and y coordinate
#' boundaries.
#'
#' @param spe The `SpatialExperiment` object representing spatially resolved transcriptomics data.
#' @param colRange a numeric vector of length two specifying the minimum and maximum x-coordinates
#' (pxl_col_in_fullres) to include in the subset. If `NULL`, the full range of x-coordinates is used.
#' @param rowRange a numeric vector of length two specifying the minimum and maximum y-coordinates
#' (pxl_row_in_fullres) to include in the subset. If `NULL`, the full range of y-coordinates is used.
#'
#' @return Returns a `SpatialExperiment` object subsetted to include only spots within the defined x and y coordinate ranges.
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
