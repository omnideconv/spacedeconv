#' Calculate Iniche recursivley
#' @param df dataframe with array_row and array_col
#' @param coordinates vector of (row, columns)
#' @param distance spot distance
#'
#' @returns list of spot ids of iniche
get_iniche <- function(df, coordinates, distance) {
  # typecheck df, coordinates and distance

  row <- coordinates[1]
  column <- coordinates[2]

  if ((row %% 2) == 0 & (column %% 2) == 1) {
    print("row and col have to be both even or odd to access spots")
    stop()
  } else {
    if ((row %% 2) == 1 & (column %% 2) == 0) {
      print("row and col have to be both even or odd to access spots")
      stop()
    }
  }

  center <- rownames(df[df$array_row == row & df$array_col == column, ])


  # if distance == 0 then return spot
  if (distance == 0) {
    return(center)
  } else {
    # else: get all surounding spots and call get_iniche for them

    # circle around the center
    s1 <- rownames(df[df$array_row == row & df$array_col == column - 2, ]) # spot below
    s2 <- rownames(df[df$array_row == row & df$array_col == column + 2, ]) # spot above
    s3 <- rownames(df[df$array_row == row - 1 & df$array_col == column + 1, ]) # spot upper right
    s4 <- rownames(df[df$array_row == row + 1 & df$array_col == column + 1, ]) # spot down right
    s5 <- rownames(df[df$array_row == row - 1 & df$array_col == column - 1, ]) # spot upper left
    s6 <- rownames(df[df$array_row == row + 1 & df$array_col == column - 1, ]) # spot upper left

    iniche <- c(s1, s2, s3, s4, s5, s6, center)

    tmp <- c(center)

    for (spot in iniche) {
      # get coordinates

      row <- df[spot, "array_row"]
      column <- df[spot, "array_col"]

      tmp <- c(tmp, get_iniche(df, c(row, column), distance - 1))
    }

    tmp <- unique(tmp)

    return(tmp)
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
  # Warning if cell types are not specified or presence and threshold are both provided
  if (is.null(cell_type_A) || is.null(cell_type_B)) {
    stop("cell type 1 or 2 missing or null")
  } else if (!(is.null(presence_matrix) & is.null(threshold))) {
    stop("only presence_marix or threshold can be supplied")
  }

  # Calculate custom presence matrix based on threshold
  if (!(is.null(threshold))) {
    presence_matrix <- presence(spe, method, threshold)
  }

  # Calculate presence matrix with an automatic method
  if (is.null(presence_matrix) & is.null(threshold)) {
    presence_matrix <- presence(spe, method)
  }

  if (distance == 0) {
    # create presence/absence vector for both celltypes
    A_pres <- presence_matrix[, cell_type_A]
    B_pres <- presence_matrix[, cell_type_B]
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

    for (i in length(iniche)) {
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
    A_pres <- niche_A_B[cell_type_A, ]
    B_pres <- niche_A_B[cell_type_B, ]
  }

  return(list(A_pres = A_pres, B_pres = B_pres))
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
#' @returns presence/absence for each cell type
#' @export

cell_pair_localization <- function(spe, method = NULL, distance = 0, density = FALSE,
                                   cell_type_A = NULL, cell_type_B = NULL,
                                   niter = 100, presence_matrix = NULL, threshold = NULL) {
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
  for (i in 1:niter) {
    coloc_rand <- vector(length = niter)
    avoid_rand <- vector(length = niter)

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
    coloc = real_coloc,
    coloc_p = p_coloc,
    coloc_rand_mean = coloc_rand_mean,
    coloc_ratio = coloc_ratio,
    avoid = real_avoid,
    avoid_p = p_avoid,
    avoid_rand_mean = avoid_rand_mean,
    avoid_ratio = avoid_ratio
  )

  # Density plot
  if (density) {
    dens <- density(coloc_rand)
    p <- plot(dens,
      main = paste0(
        "Colocalization ",
        cell_type_A, "_", cell_type_B
      ),
      xlim = range(real_coloc, coloc_rand),
      cex.axis = 1.3,
      cex.lab = 1.3,
      cex.main = 1.8
    )
    abline(v = real_coloc, col = "red")
    plot(density(avoid_rand),
      main = paste0(
        "Avoidance ",
        cell_type_A, "_", cell_type_B
      ),
      xlim = range(real_avoid, avoid_rand),
      cex.axis = 1.3,
      cex.lab = 1.3,
      cex.main = 1.8
    )
    abline(v = real_avoid, col = "red")
  }

  return(res)
}
