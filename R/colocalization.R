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


#' Calculate Colocalization
#'
#' @param spe SpatialExperiment
#' @param method presence method
#' @param distance distance in spot diameter
#' @param cell_type_1 celltype 1
#' @param cell_type_2 celltype 2
#' @param density logical
#' @param niter niter
cell_pair_localization <- function(spe, method = NULL, distance = 0, cell_type_1 = NULL, cell_type_2 = NULL, density = FALSE, niter = 100) {
  if (is.null(cell_type_1) || is.null(cell_type_2)) {
    stop("cell type 1 or 2 missing or null")
  }

  ########## cell type test

  presence <- presence(spe, method, antimode_cutoff(spe, method))
  presence[is.na(presence)] <- FALSE

  if (distance == 0) {
    # create presence/absence vector for both celltypes

    A <- presence[, cell_type_1]
    B <- presence[, cell_type_2]
  } else if (distance > 0) {
    df <- as.data.frame(SpatialExperiment::spatialCoords(spe))

    # calculate scaling factor
    scaling_offset <- 1.165 # 1.154701 # 1/cos((30/360)*2*pi) #######!!!!!!!!!!!
    smoothing_factor <- distance
    # calculate spot distance
    spot_distance <- min(sqrt((df$pxl_col_in_fullres[1] - df$pxl_col_in_fullres[-1])^2 + (df$pxl_row_in_fullres[1] - df$pxl_row_in_fullres[-1])^2)) * scaling_offset



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
    niche_pres_A <- rep(FALSE, length(iniche))
    niche_pres_B <- rep(FALSE, length(iniche))

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
    }

    # combine
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
    coloc = coloc,
    coloc_p = p_coloc,
    coloc_rand_mean = coloc_rand_mean,
    coloc_ratio = coloc_ratio,
    avoid = avoid,
    avoid_p = p_avoid,
    avoid_rand_mean = avoid_rand_mean,
    avoid_ratio = avoid_ratio
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



#' Ripley´s K celltype distribution
#'
#' @param spe SpatialExperiment
#' @param cell_type celltype of interest
#' @param method deconvolution method

ripleys_key <- function(spe, cell_type, method) {
  coords <- spatialCoords(spe)

  a <- antimode_cutoff(m = spe, method = method)
  p <- presence(spe = spe, threshold = a, method = method)
  type <- as.factor(p[, cell_type])

  pp <- spatstat.geom::ppp(
    x = coords[, 2],
    y = coords[, 1],
    xrange = range(coords[, 2]),
    yrange = range(coords[, 1]),
    marks = type
  )
  plot(spatstat.core::Kcross(pp, i = "TRUE", correction = "Ripley"), main = cell_type)
}





coloc_distance <- function(spe, method = NULL, cell_type_1, cell_type_2, distance_range = c(1, 3)){
  a <- antimode_cutoff(m = spe, method = method)
  p <- presence(spe = spe, threshold = a, method = method)
}
