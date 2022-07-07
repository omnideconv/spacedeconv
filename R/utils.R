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
checkCol <- function(object, column){
  return (column %in% names(colData(object)))
}

#' Add results to object colData
#'
#' @param spe SpatialExperiment
#' @param result deconvolution result, rows = spots, columns = cell types
addResultToObject <- function(spe, result){
  if (is.null(spe)){
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if (is.null(result)){
    stop("Parameter 'spe' is null or missing, but is required")
  }

  message("saving results to object")

  # make cell type names unique
  colnames(result) <- make.names(colnames(result))

  # check if number of spots matches result, might not be the case for all methods
  if (nrow(result) == ncol(spe)){
    # add to spatialExperiment interatively
    for (celltype in colnames(result)){
      spe[[celltype]] <- result[, celltype]
    }
  } else {
    # get the missing spots and input zero for them
    message("While saving results: dimensions don't match")

    # v1 shorter???
    # v2 longer
    v1 <- colnames(spe)
    v2 <- rownames(result)

    # get the ones from v1 missing in v2
    missing  = v1[!v1 %in% v2]

    # construct "missing data" and set all to NA
    missing_mat <- matrix(data = NA, nrow = length(missing), ncol = ncol(result))
    rownames(missing_mat) <- missing
    colnames(missing_mat) <- colnames(result)

    # construct full dataframe
    full <- rbind(result, missing_mat)

    # order accordingly
    full <- full[order(match(rownames(full), rownames(result))), ]

    # add to object
    for (celltype in colnames(full)){
      spe[[celltype]] <- full[, celltype]
    }

  }

  return (spe)

}



#' The dependencies for each method
#'
required_packages <- list(
  "autogenes" = c("reticulate"),
  "bisque" = c("BisqueRNA"),
  "bseqsc" = c("shenorrlab/bseqsc"),
  "cdseq" = c("PelzKo/CDSeq_R_Package"),
  "cibersortx" = c("uuid"),
  "cpm" = c("amitfrish/scBio"),
  "dwls" = c("PelzKo/dwls"),
  "momf" = c("grst/MOMF"),
  "music" = c("xuranw/MuSiC"),
  "scaden" = c("reticulate"),
  "scdc" = c("grst/SCDC")
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
        message(
          "To install the dependencies for all methods at once, run ",
          "devtools::install_github(\"omnideconv/omnideconv\", ",
          "dependencies = c(\"Imports\", \"Suggests\"))"
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
        message(
          "To install the dependencies for all methods at once, run ",
          "devtools::install_github(\"omnideconv/omnideconv\", ",
          "dependencies = c(\"Imports\", \"Suggests\"))"
        )
      }
      if (package_download_allowed) {
        remotes::install_github(pkgname)
      }
    }
  })
  if (repositories_set && !package_download_allowed) {
    message(
      "To install the dependencies for all methods at once, run ",
      "devtools::install_github(\"omnideconv/omnideconv\", ",
      "dependencies = c(\"Imports\", \"Suggests\"))"
    )
    stop(paste0(method, " can not be run without installing the required packages: ", packages))
  }
}
