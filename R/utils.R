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

  message("saving results to object")

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
  message("attaching token")
  if (is.null(deconvolution)) {
    stop("Deconvolution result is missing but is required")
  }

  if (is.null(token)) {
    message("No token provided, using 'deconv' as token")
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
#'
#' @export
available_results <- function(deconv) {
  if (is(deconv, "SpatialExperiment")) {
    res <- names(colData(deconv))

    res <- res[!res %in% c("in_tissue", "sample_id", "array_col", "array_row", "pxl_col_in_fullres", "pxl_row_in_fullres")]
  } else {
    print("Please provide a SpatialExperiment")
  }

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


#' Threshold a matrix
#'
#' @param spe A SpatialExperiment containing deconvolution scores for each spot and celltype
#' @param method deconvolution method
#' @param threshold if single value is provided the same threshold is used for all celltypes, it is also possible to provide a threshold vector
#'
#'
presence <- function(spe, method, threshold) {
  if (!method %in% deconvolution_methods) {
    stop("method not supported")
  }

  # create matrix with scores for each spot and celltype
  available <- available_results(spe)[startsWith(available_results(spe), method)]
  m <- as.matrix(colData(spe)[, available])

  # remove NAN????

  # calculate log(scores) +1
  m <- log(m + 1)

  # initialize matrix
  m_row <- nrow(m)
  m_col <- ncol(m)
  m_out <- matrix(FALSE, # Set all to absent (i.e. 0)
    nrow = m_row,
    ncol = m_col
  )


  # calculate threshold
  if (length(threshold) == 1) {
    m_out[m > threshold] <- TRUE
  } else if (length(threshold) == m_col) {
    for (i in 1:m_row) {
      m_out[i, ] <- m[i, ] > threshold
    }
  } else {
    stop(
      "As threshold, you can enter either a number or a vector of length ",
      m_col, "\n"
    )
  }

  rownames(m_out) <- rownames(m)
  colnames(m_out) <- colnames(m)

  return(m_out)
}

#' Determine threshold for celltype presence based on antimode of celltype density
#' @param m A SpatialExperiment Object
#' @param method deconvolution method
#' @return A vector with celltype specific cutoff values


antimode_cutoff <- function(m, method) {
  if (!method %in% deconvolution_methods) {
    stop("method not supported")
  }

  # create matrix with scores for each spot and celltype
  available <- available_results(m)[startsWith(available_results(m), method)]
  m <- colData(m)[, available]
  # threshold vector
  cutoffs <- c()
  # vector with all celltype names
  celltypes <- colnames(m)

  for (celltype in celltypes) {
    score <- m[, celltype]

    # Compute log-score

    logscore <- log(score + 1)

    # Exclude 1% most extreme values

    interval <- quantile(logscore,
      p = c(0.005, 0.995), na.rm = T
    )
    logscore <- logscore[logscore > interval[1] & logscore < interval[2]]

    # Estimate antimode

    res <- locmodes(logscore,
      mod0 = 2,
      display = F
    ) # You can put this to FALSE

    cutoff <- res$locations[2]

    # Add cutoff to threshold vector

    cutoffs <- append(cutoffs, cutoff)
  }

  names(cutoffs) <- celltypes

  return(cutoffs)
}


#' Calculate Colocalization
#'
#' @param spe SpatialExperiment
#' @param method
#' @param distance distance in spot diameter
#' @param cell_type_1
#' @param cell_type_2
#' @param density logical
#' @param niter
cell_pair_localization <- function(spe, method = NULL, distance = 0, cell_type_1 = NULL, cell_type_2 = NULL, density = FALSE, niter = 100) {
  if (is.null(cell_type_1) || is.null(cell_type_2)) {
    stop("cell type 1 or 2 missing or null")
  }

  ########## cell type test

  presence <- presence(spe, method, antimode_cutoff(spe, method))

  if (distance == 0) {
    # create presence/absence vector for both celltypes

    A <- presence[, cell_type_1]
    B <- presence[, cell_type_2]
  } else if (distance == 1) {
    df <- as.data.frame(SpatialExperiment::spatialCoords(spe))

    # calculate scaling factor
    scaling_offset <- 1.165 # 1.154701 # 1/cos((30/360)*2*pi) #######!!!!!!!!!!!
    smoothing_factor <- 1.5
    # calculate spot distance
    spot_distance <- min(sqrt((df$pxl_col_in_fullres[1] - df$pxl_col_in_fullres[-1])^2 + (df$pxl_row_in_fullres[1] - df$pxl_row_in_fullres[-1])^2)) * spot_size * scaling_offset



    # for all spots get the spots in distance and calculate mean value
    iniche <- vector(mode = "list", length = niter)
    for (spot in rownames(df)) {
      point <- df[spot, ][c("pxl_col_in_fullres", "pxl_row_in_fullres")]

      spots_in_distance <- sqrt((point$pxl_col_in_fullres - df$pxl_col_in_fullres)^2 + (point$pxl_row_in_fullres - df$pxl_row_in_fullres)^2)
      names(spots_in_distance) <- rownames(df)
      spots_in_distance <- spots_in_distance[spots_in_distance <= spot_distance * smoothing_factor]

      # new_values <- c(new_values, mean(df[names(spots_in_distance), cell_type], na.rm = TRUE))
      iniche[spot] <- list(names(spots_in_distance))
    }

    # determine presence/absence of celltypes in the iniches
    niche_pres_A <- rep(FALSE, niter)
    niche_pres_B <- rep(FALSE, niter)

    for (i in 1:length(iniche)) {
      bar <- iniche[[i]]

      uni <- presence[bar, cell_type_1]
      dui <- presence[bar, cell_type_2]

      if (sum(uni) >= 1) {
        niche_pres_A[i] <- TRUE
      }

      if (sum(dui) >= 1) {
        niche_pres_B[i] <- TRUE
      }

      niche_A_B <- rbind(niche_pres_A, niche_pres_B)
      rownames(niche_A_B) <- c(cell_type_1, cell_type_2)

      A <- niche_A_B[cell_type_1, ]
      B <- niche_A_B[cell_type_2, ]
    }

    niche_A_B <- rbind(niche_pres_A, niche_pres_B)
    rownames(niche_A_B) <- c(cell_type_1, cell_type_2)

    A <- niche_A_B[cell_type_1, ]
    B <- niche_A_B[cell_type_2, ]
  }


  loc_original <- coloc_avoid(A, B)
  coloc <- loc_original["coloc"]
  avoid <- loc_original["avoid"]

  # compute colocalization and avoidance of randomized presence matrix
  A_B_coloc_rand <- vector(length = niter)
  for (i in 1:niter) {
    A_shuffle <- sample(A)
    loc_shuffle <- coloc_avoid(A_shuffle, B)
    A_B_coloc_rand[i] <- loc_shuffle["coloc"]
  }

  A_B_avoid_rand <- vector(length = niter)
  for (i in 1:niter) {
    A_shuffle <- sample(A)
    loc_shuffle <- coloc_avoid(A_shuffle, B)
    A_B_avoid_rand[i] <- loc_shuffle["avoid"]
  }

  if (density) {
    dens <- density(A_B_coloc_rand)
    p <- plot(dens, main = paste0("Colocalization ", cell_type_1, "_", cell_type_2), xlim = range(coloc, A_B_coloc_rand))
    abline(v = coloc, col = "red")
    plot(density(A_B_avoid_rand), main = paste0("Avoidance ", cell_type_1, "_", cell_type_2), xlim = range(avoid, A_B_avoid_rand))
    abline(v = avoid, col = "red")
  }

  # p value
  p_coloc <- sum(A_B_coloc_rand > coloc) / length(A_B_coloc_rand)
  p_avoid <- sum(A_B_avoid_rand > avoid) / length(A_B_avoid_rand)

  # mean
  coloc_rand_mean <- mean(A_B_coloc_rand)
  avoid_rand_mean <- mean(A_B_avoid_rand)

  # ratio
  coloc_ratio <- coloc / coloc_rand_mean
  avoid_ratio <- avoid / avoid_rand_mean

  res <- c(
    A_B_coloc = coloc,
    A_B_coloc_p = p_coloc,
    A_B_coloc_rand_mean = coloc_rand_mean,
    A_B_coloc_ratio = coloc_ratio,
    A_B_avoid = avoid,
    A_B_avoid_p = p_avoid,
    A_B_avoid_rand_mean = avoid_rand_mean,
    A_B_avoid_ratio = avoid_ratio
  )

  return(res)
}

#' Colocalization of two celltypes
#'
#' @param A presence of cell type
#' @param B presence of cell type
coloc_avoid <- function(A, B) {
  A[is.na(A)] <- FALSE # does this make sense?????
  B[is.na(B)] <- FALSE
  coloc <- sum(A & B) / length(A)
  avoidance <- sum(!A & B) / length(A)
  return(c(coloc = coloc, avoid = avoidance))
}

#' #' @param obj SpatialExperiment
#' #' @param r rowname
#' #' @param c colname
#'
#' iniche <- function(obj, r, c){
#'
#'   if (r%%2 == 0 && c%%2 == 1){stop("Row and column have to be both even or odd")}
#'   if (c%%2 == 0 && r%%2 == 1){stop("Row and column have to be both even or odd")}
#'
#'   coord <- SpatialExperiment::spatialCoords(obj)
#'   index <- coord[coord[, "pxl_col_in_fullres"]== c, ]
#'   index <- index[index[, "pxl_row_in_fullres"] == r, ]
#' }
#'
