#' Set Path to CIBERSORT R script (`CIBERSORT.R`)
#'
#' CIBERSORT is only freely available to academic users.
#' A license an the binary can be obtained from https://cibersort.stanford.edu.
#'
#' @param path path to cibersort R script.
#'
#' @export
set_cibersort_binary <- function(path) {
  immunedeconv::set_cibersort_binary(path) # set the same for immunedeconv
  assign("cibersort_binary", path, envir = config_env)
}

#' Set Path to CIBERSORT matrix file (`LM22.txt`)
#'
#' CIBERSORT is only freely available to academic users.
#' A license an the binary can be obtained from https://cibersort.stanford.edu.
#'
#' @param path path to cibersort matrix.
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

#' Add results to object colData
#'
#' @param spatial_obj SpatialExperiment
#' @param result deconvolution result, rows = spots, columns = cell types
addResultToObject <- function(spatial_obj, result) {
  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is null or missing, but is required")
  }

  if (is.null(result)) {
    stop("Parameter 'spatial_obj' is null or missing, but is required")
  }

  # make cell type names unique
  colnames(result) <- make.names(colnames(result))

  # check if number of spots matches result, might not be the case for all methods
  if (nrow(result) == ncol(spatial_obj)) {
    # add to spatialExperiment interatively
    for (celltype in colnames(result)) {
      spatial_obj[[celltype]] <- result[, celltype]
    }
  } else {
    # get the missing spots and input zero for them
    message("While saving results: dimensions don't match")

    # v1 shorter???
    # v2 longer
    v1 <- colnames(spatial_obj)
    v2 <- rownames(result)

    # get the ones from v1 missing in v2
    missing <- v1[!v1 %in% v2]

    # construct "missing data" and set all to NA
    missing_mat <- matrix(data = NA, nrow = length(missing), ncol = ncol(result))
    rownames(missing_mat) <- missing
    colnames(missing_mat) <- colnames(result)

    # construct full dataframe
    full <- rbind(result, missing_mat)

    # order accordingly
    full <- full[order(match(rownames(full), rownames(result))), , drop = FALSE]

    # add to object
    for (celltype in colnames(full)) {
      spatial_obj[[celltype]] <- full[, celltype]
    }
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
  "cibersortx" = c("uuid"),
  "cpm" = c("amitfrish/scBio"),
  "scaden" = c("reticulate")
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
    if (bare_pkgname == "CDSeq_R_Package") {
      bare_pkgname <- "CDSeq"
    } else if (bare_pkgname == "dwls") {
      bare_pkgname <- "DWLS"
    }
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

#' Check wich deconvolutionr results are available in a SpatialExperiment object
#'
#' @param deconv SpatialExperiment
#' @param method deconvolution method (optional)
#'
#' @export
available_results <- function(deconv, method = NULL) {
  if (is(deconv, "SpatialExperiment")) {
    res <- names(colData(deconv))

    res <- res[!res %in% c("in_tissue", "sample_id", "array_col", "array_row", "pxl_col_in_fullres", "pxl_row_in_fullres")]

    if (!is.null(method)) {
      res <- res[startsWith(res, method)]
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



#' Annotate Specific Spots within a SpatialExperiment Object
#'
#' This function allows for the annotation of specified spots within a `SpatialExperiment` object.
#' Users can define a list of spots to be annotated and assign positive or negative values to these spots.
#' This is particularly useful for distinguishing or marking spots based on certain criteria or experimental results.
#'
#' @param spe A `SpatialExperiment` object representing the spatially resolved data.
#' This is the target object where the annotations will be applied.
#' @param spots A character vector or list specifying the spots within the `SpatialExperiment` object to annotate.
#' These spots should correspond to the column names of the `SpatialExperiment` object that represent specific spatial locations or features.
#' @param value_pos The value to assign to the spots specified in the `spots` parameter.
#' This value denotes a positive annotation, indicating that the spots meet a certain condition or criterion.
#' @param value_neg The value to assign to spots not specified in the `spots` parameter.
#' This value denotes a negative annotation, indicating that the spots do not meet the specified condition or criterion.
#' @param name A character string specifying the name of the annotation.
#' This name will be used to label the column in the `SpatialExperiment` object that contains the annotation values.
#'
#' @return Returns a `SpatialExperiment` object that includes the new annotations.
#' The function adds a column to the `colData` of the `SpatialExperiment` object, where the column name is specified by the `name` parameter.
#' Spots specified in the `spots` parameter are annotated with the `value_pos` value, while all other spots are annotated with the `value_neg` value.
#'
#' @export
annotate_spots <- function(spe, spots, value_pos = TRUE, value_neg = FALSE, name = "annotation") {
  df <- data.frame(row.names = colnames(spe))
  df[, name] <- value_neg
  df[spots, ] <- value_pos
  colData(spe) <- cbind(colData(spe), df)

  return(spe)
}

#' Extend a SpatialExperiment with a New Annotation Column
#'
#' This function is designed to annotate a `SpatialExperiment` object by appending a custom annotation
#' column to its `colData`. This feature allows users to integrate additional metadata or experimental
#' results that can be useful for subsequent analyses or visualizations.
#'
#' @param spatialExp A `SpatialExperiment` object that will be extended with the new annotation.
#'
#' @param columnName A `character` string that names the new annotation column to be added.
#'
#' @param values A numeric or character vector containing the annotation values to be assigned to each spot
#'
#' @return The function returns an updated `SpatialExperiment` object that includes the newly added annotation
#' column within its `colData`.
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


#' check if scaden can be found in the path variable
check_path_scaden <- function() {
  path <- paste0(reticulate::miniconda_path(), "/envs/r-omnideconv/bin/") # scaden is in there
  Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", path))
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
