#' Subset Single Cell Experiment
#'
#' Subset a SingleCellExperiment to reduce the number of cells
#'
#' @param sce SingleCellExperiment with cell type annotation
#' @param cell_type_col name of colData column containing cell type information
#' @param scenario which fraction of each celltype should be present in the subsampled dataset?
#' @param ncells number of cells in total, defaults to 1000
#' @param notEnough what to do if not enough cells are available c("remove", "asis")
#' @param seed set a seed, enables reproducibility of subsetting step
#'
#' @returns a singleCellExperiment with
#'
#' @export
#' @examples
#' data("single_cell_data_1")
#' sce <- subsetSCE(
#'   single_cell_data_1,
#'   cell_type_col = "celltype_major",
#'   scenario = "even",
#'   ncells = 100
#' )
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

  cli::cli_progress_step(msg = paste0("extracting up to ", ncells, " cells"), msg_done = paste0("Extracted ", ncells, " cells"))

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
  }


  # actually subset and return the sce
  sce <- sce[, x]

  cli::cli_progress_done()

  return(sce)
}
