#' Calculate Colocalization
#'
#' @param spe SpatialExperiment
#' @param method presence method
#' @param distance distance in spot diameter
#' @param cell_type_A celltype 1
#' @param cell_type_B celltype 2
#' @param density logical
#' @param niter permutations
#' @param presence_matrix presence matrix if available
#' @returns statistics and graph
#' @export
cell_pair_localization <- function(spe, method = NULL, distance = 0,
                                   cell_type_A = NULL, cell_type_B = NULL,
                                   density = FALSE, niter = 100, presence_matrix = NULL, threshold = NULL) {
  if (is.null(cell_type_A) || is.null(cell_type_B)) {
    stop("cell type 1 or 2 missing or null")
  }

  ########## cell type test
  # Calculate presence matrix based on abtimode cutoff approach
  if (is.null(presence_matrix)) {
    presence_matrix <- presence(spe, method)
  }
  # Set NA values to absent
  presence_matrix[is.na(presence_matrix)] <- FALSE

  # create presence/absence vector for both celltypes

  A <- presence_matrix[, cell_type_A]
  B <- presence_matrix[, cell_type_B]

  # Calculate real colocalization and avoidance events
  loc_original <- coloc_avoid(A, B)
  coloc <- loc_original["coloc"]
  avoid <- loc_original["avoid"]

  # Create randomized presence matrices and calculate corresponding colocalization/avoidance events
  coloc_rand <- vector(length = niter)
  avoid_rand <- vector(length = niter)

  for (i in 1:niter) {
    # shuffle coordinates
    presence_matrix_rand <- presence_matrix
    row.names(presence_matrix_rand) <- sample(row.names(presence_matrix_rand))

    # Create randomized presence/absence vector for both cell types
    A_rand <- presence_matrix_rand[, cell_type_A]
    B_rand <- presence_matrix_rand[, cell_type_B]

    # Calculate randomized colocalization and avoidance events
    loc_rand <- coloc_avoid(A_rand, B_rand)
    coloc_rand[i] <- loc_rand["coloc"]
    avoid_rand[i] <- loc_rand["avoid"]
  }

  if (density) {
    dens <- density(coloc_rand)
    p <- plot(dens,
      main = paste0(
        "Colocalization ",
        cell_type_A, "_", cell_type_B
      ),
      xlim = range(coloc, coloc_rand)
    )
    abline(v = coloc, col = "red")
    plot(density(avoid_rand),
      main = paste0(
        "Avoidance ",
        cell_type_A, "_", cell_type_B
      ),
      xlim = range(avoid, avoid_rand)
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





##### Add distance parameter
##
if (distance == 0) {

} else if (distance > 0) {
  df <- as.data.frame(colData(spe))

  # # for all spots get the spots in distance and calculate mean value
  iniche <- vector(mode = "list", length = niter)
  for (spot in rownames(df)) {
    iniche[spot] <- list(get_iniche(df, get_spot_coordinates(df, spot), distance = distance))
  }

  # determine presence/absence of celltypes in the iniches
  niche_pres_A <- rep(FALSE, length(iniche))
  niche_pres_B <- rep(FALSE, length(iniche))

  for (i in 1:length(iniche)) {
    bar <- iniche[[i]]

    uni <- presence[bar, cell_type_A]
    dui <- presence[bar, cell_type_B]

    if (sum(uni) >= 1) {
      niche_pres_A[i] <- TRUE
    }

    if (sum(dui) >= 1) {
      niche_pres_B[i] <- TRUE
    }
  }

  # combine
  niche_A_B <- rbind(niche_pres_A, niche_pres_B)
  rownames(niche_A_B) <- c(cell_type_A, cell_type_B)

  A <- niche_A_B[cell_type_A, ]
  B <- niche_A_B[cell_type_B, ]
}
##
