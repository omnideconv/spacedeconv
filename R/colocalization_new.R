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

  # Warning if cell types are not specified or presence and threshold are both provided
  if (is.null(cell_type_A) || is.null(cell_type_B)) {
    stop("cell type 1 or 2 missing or null")
  }
   else if (!(is.null(presence_matrix) & is.null(threshold))) {
    stop("only presence_marix or threshold can be supplied")
  }

  # Calculate custom presence matrix based on threshold
  if (!(is.null(threshold))) {
    presence_matrix <- presence(spe, method, threshold)
  }

  # Calculate presence matrix based with an automatic method
  if (is.null(presence_matrix) & is.null(threshold)) {
    presence_matrix <- presence(spe, method)
  }

  if (distance == 0) {

    # create presence/absence vector for both celltypes
    cellA_pres <- presence_matrix[, cell_type_A]
    cellB_pres <- presence_matrix[, cell_type_B]

  } else if (distance > 0) {
    df <- as.data.frame(colData(spe))

    # for all spots get the spots in distance and calculate mean value
    iniche <- vector(mode = "list", length = niter)
    for (spot in rownames(df)) {
      iniche[spot] <- list(get_iniche(df, get_spot_coordinates(df, spot), distance = distance))
    }

    # determine presence/absence of celltypes in the iniches
    niche_pres_A <- rep(FALSE, length(iniche))
    niche_pres_B <- rep(FALSE, length(iniche))

    for (i in 1:length(iniche)) {
      # extract spots of iniche
      bar <- iniche[[i]]
      # calculate presence values for cell type A and B in the iniche
      uni <- presence[bar, cell_type_A]
      dui <- presence[bar, cell_type_B]
      # set whole iniche value to present, if at least one spot contains cell type A
      if (sum(uni) >= 1) {
        niche_pres_A[i] <- TRUE
      }
      # set whole iniche value to present, if at least one spot contains cell type B
      if (sum(dui) >= 1) {
        niche_pres_B[i] <- TRUE
      }
    }

    # combine presence/absence values for both cell types
    niche_A_B <- rbind(niche_pres_A, niche_pres_B)
    rownames(niche_A_B) <- c(cell_type_A, cell_type_B)

    # create presence/absence vectors for cell type A and B
    cellA_pres <- niche_A_B[cell_type_A, ]
    cellB_pres <- niche_A_B[cell_type_B, ]
  }

  return(c(cellA_pres = cellA_pres, cellB_pres = cellB_pres))
}


#' Colocalization of two celltypes
#'
#' @param A presence of cell type
#' @param B presence of cell type
#'
coloc_avoid <- function(A, B) {
  coloc <- sum(A & B) / length(A)
  avoidance <- sum(!A & B) / length(A)
  return(c(coloc = coloc, avoid = avoidance))
}



  # Calculate real colocalization and avoidance events based on the presence/absence vectors for cell type A and B
  loc_original <- coloc_avoid(cellA_pres, cellB_pres)
  coloc <- loc_original["coloc"]
  avoid <- loc_original["avoid"]



  #################  Randomization??
  # Create randomized presence matrices and calculate corresponding colocalization/avoidance events
  coloc_rand <- vector(length = niter)
  avoid_rand <- vector(length = niter)

  for (i in 1:niter) {
    # shuffle coordinates of presence/absence vectors --> names need to be in the same order for coloc_avoid function
    A_rand <- sample(cellA_pres)
    names(A_rand) <- names(cellA_pres)
    B_rand <- sample(cellB_pres)
    names(A_rand) <- names(cellB_pres)

    # Calculate presence matrix for randomized version
    presence_matrix_rand <- presence(spe, method)


    # Calculate randomized colocalization and avoidance events
    loc_rand <- coloc_avoid(A_rand, B_rand)
    coloc_rand[i] <- loc_rand["coloc"]
    avoid_rand[i] <- loc_rand["avoid"]
  }

  #################  Increase title size and size of axis values!
  if (density) {
    dens <- density(coloc_rand)
    p <- plot(dens,
      main = paste0(
        "Colocalization ",
        cell_type_A, "_", cell_type_B
      ),
      xlim = range(coloc, coloc_rand),
      cex.axis = 1.3,
      cex.lab = 1.3,
      cex.main = 1.8
    )
    abline(v = coloc, col = "red")
    plot(density(avoid_rand),
      main = paste0(
        "Avoidance ",
        cell_type_A, "_", cell_type_B
      ),
      xlim = range(avoid, avoid_rand),
      cex.axis = 1.3,
      cex.lab = 1.3,
      cex.main = 1.8
    )
    abline(v = avoid, col = "red")
  }

  # p value
  p_coloc <- sum(coloc_rand > coloc) / length(coloc_rand)
  p_avoid <- sum(avoid_rand > avoid) / length(avoid_rand)

  # mean
  coloc_rand_mean <- mean(coloc_rand)
  avoid_rand_mean <- mean(avoid_rand)

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
