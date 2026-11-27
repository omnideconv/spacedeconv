#' Register the Path to the CIBERSORT Script
#'
#' Stores the path to the `CIBERSORT.R` script in the spacedeconv configuration and
#' forwards it to `immunedeconv::set_cibersort_binary()` so CIBERSORT-based methods
#' can run from this R session.
#'
#' CIBERSORT is only freely available to academic users. A license and script can
#' be obtained from https://cibersort.stanford.edu.
#'
#' @param path Path to the `CIBERSORT.R` script.
#'
#' @export
set_cibersort_binary <- function(path) {
  immunedeconv::set_cibersort_binary(path) # set the same for immunedeconv
  assign("cibersort_binary", path, envir = config_env)
}

#' Register the Path to the CIBERSORT Signature Matrix
#'
#' Stores the path to the `LM22.txt` signature matrix in the spacedeconv configuration
#' and forwards it to `immunedeconv::set_cibersort_mat()` so CIBERSORT-based methods
#' can access the signature from this R session.
#'
#' CIBERSORT is only freely available to academic users. A license and matrix can
#' be obtained from https://cibersort.stanford.edu.
#'
#' @param path Path to the `LM22.txt` signature matrix file.
#'
#' @export
set_cibersort_mat <- function(path) {
  immunedeconv::set_cibersort_mat(path)
  assign("cibersort_mat", path, envir = config_env)
}

#' get matrices from SingleCellExperiment
#' @param single_cell_object SingleCellExperiment
#' @param cell_type_col column containing the cell type
#' @returns list containing the expression matrix and cell type annotation vector

getMatricesFromSCE <- function(single_cell_object, cell_type_col = "cell_ontology_class") {
  # tests

  # test if object valid

  # test if cell_type_col in names(colData())

  counts <- as(counts(single_cell_object), "dgCMatrix") # count matrix as sparse matrix
  cell_type_annotation <- as.character(colData(single_cell_object)[[cell_type_col]])


  return(list(counts = counts, cell_type_annotation = cell_type_annotation))
}

#' Check if Column exists in object
#' @param object SingleCellExperiment or SpatialExperiment
#' @param column column name to check for existence
#' @returns if column exists in object
checkCol <- function(object, column) {
  return(column %in% names(SingleCellExperiment::colData(object)))
}

#' Import External Deconvolution Results
#'
#' Adds externally computed cell-type estimates to an existing SpatialExperiment.
#' Rows are aligned by spot ID. Values are preserved without normalization;
#' expression data, spatial coordinates, and images are unchanged.
#'
#' @param spe A `SpatialExperiment` with unique, non-empty spot IDs.
#' @param results A numeric matrix or data frame with spots as rows and cell types
#' as columns. Row names must match all `colnames(spe)` exactly, in any order.
#' Values must be finite and non-negative; they need not sum to one.
#' @param result_name A prefix identifying the external result, for example
#' `"omnideconv_bayesprism"`. Use a syntactically valid R name starting with a
#' letter. The names `"expression"` and `"cluster"` are reserved.
#' @return The updated `SpatialExperiment`. Results are stored in `colData` as
#' `<result_name>_<cell_type>`, with cell-type names converted using `make.names()`.
#' Ambiguous names and collisions with existing columns are rejected.
#' Use `result_name` with `available_results()`, `plot_spatial()`,
#' and `cluster()` to select the imported result group.
#' @export
import_deconvolution_results <- function(spe, results, result_name) {
  if (!is(spe, "SpatialExperiment")) {
    stop("spe must be a SpatialExperiment.", call. = FALSE)
  }
  if (!is.character(result_name) || length(result_name) != 1L || is.na(result_name) ||
      !grepl("^[A-Za-z][A-Za-z0-9._]*$", result_name) ||
      make.names(result_name) != result_name || result_name %in% c("expression", "cluster")) {
    stop("result_name must be a valid R name starting with a letter, other than 'expression' or 'cluster'.", call. = FALSE)
  }
  if (!(is.matrix(results) || is.data.frame(results)) ||
      nrow(results) == 0L || ncol(results) == 0L) {
    stop("results must be a non-empty matrix or data frame (rows = spots, columns = cell types).", call. = FALSE)
  }
  results <- as.matrix(results)
  if (!is.numeric(results) || any(!is.finite(results)) || any(results < 0)) {
    stop("results must contain finite, non-negative numeric values.", call. = FALSE)
  }
  valid_names <- function(x) {
    !is.null(x) && !anyNA(x) && all(nzchar(trimws(x))) && !anyDuplicated(x)
  }
  if (!valid_names(colnames(spe)) || !valid_names(rownames(results)) ||
      !setequal(rownames(results), colnames(spe))) {
    stop("Result row names must match all unique, non-empty spot IDs in spe exactly.", call. = FALSE)
  }
  cell_types <- colnames(results)
  if (!valid_names(cell_types) || anyDuplicated(make.names(cell_types))) {
    stop("Cell-type column names must be non-empty and unique, including after make.names() conversion.", call. = FALSE)
  }
  colnames(results) <- paste0(result_name, "_", make.names(cell_types))
  collisions <- intersect(colnames(results), names(colData(spe)))
  if (length(collisions) > 0L) {
    stop(paste("Result columns already exist:", paste(collisions, collapse = ", ")),
      call. = FALSE)
  }
  addResultToObject(spe, results)
}

#' Add results to object colData
#'
#' @param spatial_obj SpatialExperiment
#' @param result deconvolution result, rows = spots, columns = cell types
addResultToObject <- function(spatial_obj, result) {
  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is null or missing, but is required")
  }

  if (is.null(result)) {
    stop("Parameter 'result' is null or missing, but is required")
  }

  result <- as.matrix(result)
  ids <- rownames(result)
  if (is.null(ids) || anyNA(ids) || anyDuplicated(ids)) {
    stop("Result spot IDs must be present and unique")
  }
  if (any(!ids %in% colnames(spatial_obj))) {
    stop("Result contains unknown spot IDs")
  }
  colnames(result) <- make.names(colnames(result), unique = TRUE)
  # Match even when dimensions agree: backends may reorder or omit spots.
  aligned <- result[match(colnames(spatial_obj), ids), , drop = FALSE]
  for (celltype in colnames(aligned)) {
    spatial_obj[[celltype]] <- unname(aligned[, celltype])
  }

  return(spatial_obj)
}

#' get deconvolution results from object
#' @param spatial_obj SpatialExperiment
get_results_from_object <- function(spatial_obj) {
  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is null or missing, but is required")
  }

  tmp <- SingleCellExperiment::colData(spatial_obj)
  tmp <- as.matrix(tmp[, -1]) # ???

  return(tmp)
}

#' The dependencies for each method
#'
required_packages <- list(
  "cpm" = c("amitfrish/scBio")
)


#' Checking and installing all dependencies for the specific methods
#'
#' @param method The name of the method that is used
check_and_install <- function(method) {
  if (!(method %in% deconvolution_methods)[[1]]) {
    stop(
      paste(
        "Method", method,
        "not recognized. Please refer to 'deconvolution_methods' for the integrated methods."
      )
    )
  }
  method <- method[[1]]
  packages <- required_packages[[method]]
  github_pkgs <- grep("^.*?/.*?$", packages, value = TRUE)
  cran_pkgs <- packages[!(packages %in% github_pkgs)]
  repositories_set <- FALSE
  package_download_allowed <- FALSE
  sapply(cran_pkgs, function(pkgname) {
    if (!requireNamespace(pkgname, quietly = TRUE)) {
      if (!repositories_set) {
        utils::setRepositories(graphics = FALSE, ind = c(1, 2, 3, 4, 5))
        repositories_set <<- TRUE
        package_download_allowed <<- askYesNo(
          paste0(
            "You requested to run ", method,
            " which is currently not installed. Do you want ",
            "to install the packages required for it: ", packages
          )
        )
      }
      if (package_download_allowed) {
        utils::install.packages(pkgname)
      }
    }
  })
  sapply(github_pkgs, function(pkgname) {
    bare_pkgname <- sub(".*?/", "", pkgname)
    if (!requireNamespace(bare_pkgname, quietly = TRUE)) {
      if (!repositories_set) {
        utils::setRepositories(graphics = FALSE, ind = c(1, 2, 3, 4, 5))
        repositories_set <<- TRUE
        package_download_allowed <<- askYesNo(
          paste0(
            "You requested to run ", method,
            " which is currently not installed. Do you want ",
            "to install the packages required for it: ", packages
          )
        )
      }
      if (package_download_allowed) {
        remotes::install_github(pkgname)
      }
    }
  })
  if (repositories_set && !package_download_allowed) {
    stop(paste0(method, " can not be run without installing the required packages: ", packages))
  }
}

#' Remove Spots with zero expression
#'
#' This function removes spots/columns with zero expression detected. These spots might result in errors during computation
#'
#' @param object SummarizedExperiment or any related datatypes
#'
#' @returns Expression Object without all zero columns
removeZeroExpression <- function(object) {
  # ensure that library size > 0
  nspots <- sum(Matrix::colSums(counts(object)) == 0)
  if (nspots > 0) {
    # remove spots with all zero expression
    message("removing ", nspots, " spots with zero expression")
    object <- object[, !Matrix::colSums(counts(object)) == 0]
  }

  return(object)
}

#' Check Rowname/Colname Presence
#'
#' Check for Rowname and Column Name existence in expression objects
#'
#' @param object SingleCellExperiment or SpatialExperiment
#'
#' @returns boolean, TRUE if one of rownames/colnames is NULL
checkRowColumn <- function(object) {
  return(is.null(rownames(object)) || is.null(colnames(object)))
}

#' Attach method token to deconvolution result
#'
#' Rename Celltypes of deconvolution result and add method token
#' @param deconvolution deconvolution result as matrix
#' @param token method name or custom token
#'
#' @returns deconvolution result with renamed celltypes
attachToken <- function(deconvolution, token = "deconv") {
  if (is.null(deconvolution)) {
    stop("Deconvolution result is missing but is required")
  }

  # get colnames, attach token and overwrite
  celltypes <- colnames(deconvolution)
  celltypes <- paste0(token, "_", celltypes)
  colnames(deconvolution) <- celltypes

  return(deconvolution)
}

#' List Available Deconvolution Results in a SpatialExperiment
#'
#' Returns column names in `colData` that correspond to deconvolution results.
#' Optionally filters by a method prefix (e.g., "spatialdwls").
#'
#' @param deconv A `SpatialExperiment` containing deconvolution results.
#' @param method Optional prefix used to filter result columns (typically the
#' internal method token from `spacedeconv::deconvolution_methods`, or a custom
#' `result_name` used when running or importing a method). Prefixes are matched
#' at the underscore separating the method name from the result name.
#'
#' @export
available_results <- function(deconv, method = NULL) {
  if (is(deconv, "SpatialExperiment")) {
    res <- names(colData(deconv))

    res <- res[!res %in% c("in_tissue", "sample_id", "array_col", "array_row", "pxl_col_in_fullres", "pxl_row_in_fullres")]

    if (!is.null(method)) {
      res <- res[startsWith(res, paste0(sub("_$", "", method), "_"))]
    }
  } else {
    print("Please provide a SpatialExperiment")
  }

  res <- sort(res)

  return(res)
}

#' Check for ENSEBL IDs
#'
#' @param names vector of rownames
#'
#' @returns TRUE if all are ensembl
checkENSEMBL <- function(names) {
  if (sum(grepl("^ENS", names)) / length(names) >= 0.05) {
    # more than 5% are ENSEMBL
    return(TRUE)
  } else {
    return(FALSE)
  }
}



#' Annotate Specific Spots in a SpatialExperiment
#'
#' Adds a new `colData` column that marks selected spots with `value_pos` and all
#' other spots with `value_neg`.
#'
#' @param spe A `SpatialExperiment` to annotate.
#' @param spots Character vector of spot IDs (column names) to mark.
#' @param value_pos Value assigned to selected spots.
#' @param value_neg Value assigned to all other spots.
#' @param name Name of the new annotation column.
#'
#' @return The updated `SpatialExperiment` with the new annotation column.
#'
#' @export
annotate_spots <- function(spe, spots, value_pos = TRUE, value_neg = FALSE, name = "annotation") {
  df <- data.frame(row.names = colnames(spe))
  df[, name] <- value_neg
  df[spots, ] <- value_pos
  colData(spe) <- cbind(colData(spe), df)

  return(spe)
}

#' Add a Custom Annotation Column
#'
#' Appends a new column to `colData` with user-provided values. This can be
#' used to add metadata for downstream analyses or visualization.
#'
#' @param spatialExp A `SpatialExperiment` to annotate.
#' @param columnName Name of the new annotation column.
#' @param values Vector of annotation values (length must match `ncol(spatialExp)`).
#'
#' @return The updated `SpatialExperiment` with the new annotation column.
#'
#' @export
addCustomAnnotation <- function(spatialExp, columnName, values) {
  # Checking if the length of values matches the number of columns in spatialExp
  if (length(values) != ncol(spatialExp)) {
    stop("The length of values must be equal to the number of columns in the SpatialExperiment.")
  }

  # Adding the new column to colData
  colData(spatialExp)[[columnName]] <- values

  # Return the updated SpatialExperiment
  return(spatialExp)
}

#' Get coordinates of spot id
#' @param df colData dataframe
#' @param spotid spotid
get_spot_coordinates <- function(df, spotid) {
  df <- as.data.frame(df)
  return(c(df[spotid, "array_row"], df[spotid, "array_col"]))
}

#' Convert assays to sparse Matrices
#'
#' @param spe SpatialExperiment
#' @param assay assay to use
check_datatype <- function(spe, assay = "counts") {
  if (!is(assay(spe, assay), "dgCMatrix") && !is(assay(spe, assay), "TENxMatrix")) {
    assay(spe, assay) <- as(assay(spe, assay), "sparseMatrix")

    cli::cli_alert_info("Converting data to sparse matrices")
  }

  return(spe)
}
