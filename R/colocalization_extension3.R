#' Heatmaps cell pair localization and correlation of scores
#'
#' Generates a heatmap of cell pair localization ratios and p-values as asterisks as well as a correlation heatmap of scores and correlation probabilities as asterisks between all cell types found after deconvolution
#'
#' @param spe SpatialExperiment object
#' @param method deconvolution method
#' @param distance size area of interest
#' @param correlation correlation heatmap based on scores
#' @param localization colocalization and avoidance heatmap
#' @param matrix return colocalization, avoidance, correlation matrices separately
#'
#' @returns cell pair localization heatmap and correlation heatmap
#'
#' @export

localization_heatmap <- function(spe, method, distance = 0, correlation = TRUE, localization = TRUE, matrix = FALSE) {
  # create matrix with scores for each spot and celltype
  available <- available_results(spe)[startsWith(available_results(spe), method)]
  m <- as.matrix(colData(spe)[, available])

  if (localization == TRUE){
    # Create empty matrix to fill in p-values
    p_coloc <- matrix(ncol = ncol(m), nrow = ncol(m))
    p_avoid <- matrix(ncol = ncol(m), nrow = ncol(m))
    rownames(p_coloc) <- colnames(m)
    colnames(p_coloc) <- colnames(m)
    rownames(p_avoid) <- colnames(m)
    colnames(p_avoid) <- colnames(m)
    ratio_coloc <- matrix(ncol = ncol(m), nrow = ncol(m))
    ratio_avoid <- matrix(ncol = ncol(m), nrow = ncol(m))
    colnames(ratio_coloc) <- colnames(m)
    rownames(ratio_coloc) <- colnames(m)
    colnames(ratio_avoid) <-colnames(m)
    rownames(ratio_avoid) <- colnames(m)


    presence <- presence(spe, method)

    # Calculate pairwise p-values and ratios
    for (i in 1:length(rownames(p_coloc))) {
      for (j in 1:length(colnames(p_coloc))) {
        loco <- cell_pair_localization(spe, method = method, cell_type_1 = rownames(p_coloc)[i], cell_type_2 = rownames(p_coloc)[j], density = F, distance = distance, presence = presence)
        p_coloc[i, j] <- loco["coloc_p"]
        p_avoid[i, j] <- loco["avoid_p"]
        ratio_coloc[i,j] <- loco["coloc_ratio.coloc"]
        ratio_avoid[i,j] <- loco["avoid_ratio.avoid"]
      }
    }

    # Colocalisation ration with p-value
    corrplot(ratio_coloc, p.mat = p_coloc, method = "color", diag = FALSE, type = "lower",
             is.corr = FALSE,
             col = COL1("Blues"),
             sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
             insig = "label_sig", pch.col = "black", tl.col = "black",
             title = "Colocalization ratio with significance",
             mar = c(0,0,2,0))

    # Avoidance ratio with p-value
    corrplot(ratio_avoid, p.mat = p_avoid, method = "color", diag = FALSE, type = "lower",
           is.corr = FALSE,
           col = COL1("Reds"),
           sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
           insig = "label_sig", pch.col = "black", tl.col = "black",
           title = "Avoidance ratio with significance",
           mar = c(0,0,2,0))
  cat(paste0("\t", c("Significance Colocalization/Avoidance:", "* <0.05", "** <0.01", "*** <0.001"), "\n"))

  if (matrix == TRUE) {
    print("Ratio_coloc")
    print(ratio_coloc)
    print("p_coloc")
    print(p_coloc)
    print("Ratio_avoid")
    print(ratio_avoid)
    print("p_avoid")
    print(p_avoid)
  }
  }

  if (correlation == TRUE){
    correlation <- corr.test(m, adjust = "none")$r
    cor_prob <- corr.test(m, adjust = "none")$p
    corrplot(correlation, p.mat = cor_prob, method = 'color', diag = FALSE, type = 'lower',
             sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
             insig = 'label_sig', pch.col = 'black', tl.col = "black",
             title = "Correlation with significance",
             mar=c(0,0,2,0))
    cat(paste0("\t", c("Significance correlation:", "* <0.05", "** <0.01", "*** <0.001"), "\n"))

    if (matrix == TRUE) {
      print("Correlation of scores")
      print(correlation)
      print("p_correlation")
      print(cor_prob)
    }
  }
}


