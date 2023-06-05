#' Threshold a matrix
#'
#' @param spe A SpatialExperiment containing deconvolution scores for each spot and celltype
#' @param method deconvolution method
#' @param threshold if single value is provided the same threshold is used for all celltypes, it is also possible to provide a threshold vector
#'
#' @export
presence <- function(spe, method, threshold = NULL) {
  if (!method %in% deconvolution_methods) {
    stop("method not supported")
  }

  # create matrix with scores for each spot and celltype
  available <- available_results(spe)[startsWith(available_results(spe), method)]
  m <- as.matrix(colData(spe)[, available])

  # remove NAN????
  m[is.nan(m)] <- 0

  # calculate log(scores) +1
  m <- log(m + 1)

  # initialize matrix
  m_row <- nrow(m)
  m_col <- ncol(m)
  m_out <- matrix(FALSE, # Set all to absent (i.e. 0)
    nrow = m_row,
    ncol = m_col
  )


  if (is.null(threshold)) {
    print("Calculating Antimode cutoffs")
    threshold <- antimode_cutoff(spe, method)
  }


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
#' @param spe A SpatialExperiment Object
#' @param method deconvolution method
#' @return A vector with celltype specific cutoff values


antimode_cutoff <- function(spe, method) {
  if (!method %in% deconvolution_methods) {
    stop("method not supported")
  }

  # create matrix with scores for each spot and celltype
  available <- available_results(spe)[startsWith(available_results(spe), method)]
  spe <- colData(spe)[, available]
  # threshold vector
  cutoffs <- c()
  # vector with all celltype names
  celltypes <- colnames(spe)

  for (celltype in celltypes) {
    score <- spe[, celltype]

    # Compute log-score

    logscore <- log(score + 1)

    # Exclude 1% most extreme values

    interval <- quantile(logscore,
      p = c(0.005, 0.995), na.rm = T
    )

    ###### PARAMETER ohne PROZENT, default hardcoded 1%
    logscore <- logscore[logscore > interval[1] & logscore < interval[2]]
    ###### PARAMETER

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



#' Calculate Iniche recursivley
#' @param df dataframe with array_row and array_col
#' @param coordinates vector of (row, columns)
#' @param distance spot distance
#'
#' @returns list of spot ids of iniche

get_iniche <- function(df, coordinates, distance) {
  # extract coordinates from the input
  row <- coordinates[1]
  column <- coordinates[2]

  # check format of coordinate input
  if ((row %% 2) == 0 & (column %% 2) == 1) {
    print("row and col have to be both even or odd to access spots")
    stop()
  } else {
    if ((row %% 2) == 1 & (column %% 2) == 0) {
      print("row and col have to be both even or odd to access spots")
      stop()
    }
  }

  # Extract the center spot based on given coodinates
  center <- rownames(df[df$array_row == row & df$array_col == column, ])


  # if distance == 0 then return spot
  if (distance == 0) {
    return(center)
  } else {
    # else: get all surounding spots of the center spot and call get_iniche for them

    # circle around the center to get the iniche of the center spot
    s1 <- rownames(df[df$array_row == row & df$array_col == column - 2, ]) # spot below
    s2 <- rownames(df[df$array_row == row & df$array_col == column + 2, ]) # spot above
    s3 <- rownames(df[df$array_row == row - 1 & df$array_col == column + 1, ]) # spot upper right
    s4 <- rownames(df[df$array_row == row + 1 & df$array_col == column + 1, ]) # spot down right
    s5 <- rownames(df[df$array_row == row - 1 & df$array_col == column - 1, ]) # spot upper left
    s6 <- rownames(df[df$array_row == row + 1 & df$array_col == column - 1, ]) # spot upper left

    iniche <- c(s1, s2, s3, s4, s5, s6, center)

    spots_group <- c(center)

    # Get iniche for every spot in the iniche of the center spot
    for (spot in iniche) {
      # get coordinates
      row <- df[spot, "array_row"]
      column <- df[spot, "array_col"]

      spots_group <- c(spots_group, get_iniche(df, c(row, column), distance - 1))
    }

    spots_group <- unique(spots_group)

    return(spots_group)
  }
}

#' Calculate Colocalization
#'
#' @param spe SpatialExperiment
#' @param method presence method
#' @param distance distance in spot diameter
#' @param cell_type_A celltype 1
#' @param cell_type_B celltype 2
#' @param niter permutations
#' @param presence_matrix presence matrix if available
#' @param threshold threshold for custom presence determination
#' @returns presence/absence for each cell type
#' @export

cell_pair_presence <- function(spe, method = NULL, distance = 0,
                               cell_type_A = NULL, cell_type_B = NULL,
                               niter = 100, presence_matrix = NULL, threshold = NULL) {
  # Stop if cell types are not specified
  if (is.null(cell_type_A) || is.null(cell_type_B)) {
    stop("'cell_type_A' type 'cell_type_B' must be provided")
  }

  # IF presence_matrix is provided, do nothing
  # IF not provided, presence must be computed
  if (is.null(presence_matrix)) {
    if (is.null(method)) {
      stop("If 'presence_matrix' is not provided, 'method' needs to be specified")
    } else {
      message("Calculating presence matrix...")
      presence_matrix <- presence(spe, method, threshold)
    }
  } else {
    message("Using provided presence matrix...")
  }

  if (distance < 0) {
    stop("'distance' must be non-negative")
  } else if (distance == 0) {
    # create presence/absence vector for both celltypes
    A_pres <- presence_matrix[, cell_type_A]
    B_pres <- presence_matrix[, cell_type_B]
  } else {
    # FF: Suggestion, don't use `presence()` again (especially because you are
    # not even passing the parameters), but compute the niche coordinates
    # and use them to derive niche presence from `presence_matrix` computed above

    df <- as.data.frame(colData(spe))

    # for all spots get the spots in distance and calculate mean value
    iniche <- vector(mode = "list", length = length(rownames(df)))
    for (i in 1:length(iniche)) {
      iniche[[i]] <- get_iniche(df, get_spot_coordinates(df, rownames(df)[i]), distance = distance)
    }

    # determine presence/absence of celltypes in the iniches
    A_pres <- rep(FALSE, length(iniche))
    B_pres <- rep(FALSE, length(iniche))

    for (i in 1:length(iniche)) {
      # extract presence values for iniche
      pres_niche <- presence_matrix[iniche[[i]], ]
      # set whole iniche value to present, if at least one spot contains cell type A
      if (sum(pres_niche[, cell_type_A]) >= 1) {
        A_pres[i] <- TRUE
      }
      # set whole iniche value to present, if at least one spot contains cell type B
      if (sum(pres_niche[, cell_type_B]) >= 1) {
        B_pres[i] <- TRUE
      }
    }
  }

  return(list(A_pres = A_pres, B_pres = B_pres))
}


#' Colocalization of two celltypes
#'
#' @param A presence of cell type
#' @param B presence of cell type
#'
coloc_avoid <- function(A, B) {
  ABsum <- A + B
  coloc <- sum(ABsum == 2) / length(ABsum)
  avoid <- sum(ABsum == 1) / length(ABsum)
  return(c(coloc = coloc, avoid = avoid))
}

#' Cell pair colocalization or avoidance
#' @param spe SpatialExperiment
#' @param method presence method
#' @param density density plot
#' @param distance distance in spot diameter
#' @param cell_type_A celltype 1
#' @param cell_type_B celltype 2
#' @param niter permutations
#' @param presence_matrix presence matrix if available
#' @param threshold threshold for custom presence determination
#' @param title adjust title of density plot
#' @returns presence/absence for each cell type
#' @export

cell_pair_localization <- function(spe, method = NULL, distance = 0, density = FALSE,
                                   cell_type_A = NULL, cell_type_B = NULL,
                                   niter = 100, presence_matrix = NULL, threshold = NULL, title = NULL) {
  # Calculate cell type presence
  pair_pres <- cell_pair_presence(spe,
    method = method, distance = distance, cell_type_A = cell_type_A, cell_type_B = cell_type_B,
    niter = niter, presence_matrix = presence_matrix, threshold = threshold
  )
  cellA_pres <- pair_pres$A_pres
  cellB_pres <- pair_pres$B_pres

  # Calculate real colocalization and avoidance events
  real_coloc <- coloc_avoid(cellA_pres, cellB_pres)["coloc"]
  real_avoid <- coloc_avoid(cellA_pres, cellB_pres)["avoid"]

  # Randomize presence/absence vectors independently and determine colocalization/avoidance events
  coloc_rand <- vector(length = niter)
  avoid_rand <- vector(length = niter)
  for (i in 1:niter) {
    # shuffle coordinates of presence/absence vectors --> names need to be in the same order for coloc_avoid function
    A_rand <- sample(cellA_pres)
    names(A_rand) <- names(cellA_pres)
    B_rand <- sample(cellB_pres)
    names(B_rand) <- names(cellB_pres)

    # Calculate colocalization/avoidance events of randomized version
    coloc_rand[i] <- coloc_avoid(A_rand, B_rand)["coloc"]
    avoid_rand[i] <- coloc_avoid(A_rand, B_rand)["avoid"]
  }

  # Colocalization/avoidance statistics

  # p value
  p_coloc <- sum(coloc_rand > real_coloc) / length(coloc_rand)
  p_avoid <- sum(avoid_rand > real_avoid) / length(avoid_rand)

  # mean
  coloc_rand_mean <- mean(coloc_rand)
  avoid_rand_mean <- mean(avoid_rand)

  # ratio
  coloc_ratio <- real_coloc / coloc_rand_mean
  avoid_ratio <- real_avoid / avoid_rand_mean

  res <- c(
    coloc_real = round(real_coloc, digits = 2),
    coloc_p = round(p_coloc, digits = 2),
    coloc_rand_mean = round(coloc_rand_mean, digits = 2),
    coloc_ratio = round(coloc_ratio, digits = 2),
    avoid_real = round(real_avoid, digits = 2),
    avoid_p = round(p_avoid, digits = 2),
    avoid_rand_mean = round(avoid_rand_mean, digits = 2),
    avoid_ratio = round(avoid_ratio, digits = 2)
  )

  # Density plot
  if (density) {
    # Title
    if (is.null(title)) {
      title_coloc <- paste0(
        "Colocalization ",
        cell_type_A, "_", cell_type_B
      )
      title_avoid <- paste0(
        "Avoidance ",
        cell_type_A, "_", cell_type_B
      )
    } else if (!(is.null(title))) {
      title_coloc <- paste0("Colocalization", "_", title)
      title_avoid <- paste0("Avoidance", "_", title)
    }
    # Colocalization plot
    dens <- density(coloc_rand)
    p <- plot(dens,
      main = title_coloc,
      xlim = range(real_coloc, coloc_rand),
      cex.axis = 1.3,
      cex.lab = 1.3,
      cex.main = 1.8
    )
    abline(v = real_coloc, col = "red")

    # Avoidance plot
    plot(density(avoid_rand),
      main = title_avoid,
      xlim = range(real_avoid, avoid_rand),
      cex.axis = 1.3,
      cex.lab = 1.3,
      cex.main = 1.8
    )
    abline(v = real_avoid, col = "red")
  }

  return(res)
}
