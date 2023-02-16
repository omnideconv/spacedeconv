#' Heatmaps cell pair localization and correlation of scores
#'
#' Generates a heatmap of cell pair localization p-values and a correlation heatmap of scores between all celltypes found after deconvolution
#'
#' @param spe SpatialExperiment object
#' @param method deconvolution method
#' @param distance size area of interest
#' @param correlation correlation heatmap based on scores
#' @param localization colocalization and avoidance heatmap
#' @param matrix return colocalization, avoidance, correlation matrices separately
#' @param alpha significance niveau
#'
#' @returns cell pair localization heatmap and correlation heatmap
#'
#' @export

localization_heatmap <- function(spe, method, distance = 0, correlation = TRUE, localization = TRUE, matrix = FALSE, alpha = 0.05) {
  # create matrix with scores for each spot and celltype
  available <- available_results(spe)[startsWith(available_results(spe), method)]
  m <- as.matrix(colData(spe)[, available])

  if (localization == TRUE) {
    # Create empty matrix to fill in p-values
    mat_coloc <- matrix(ncol = ncol(m), nrow = ncol(m))
    mat_avoid <- matrix(ncol = ncol(m), nrow = ncol(m))
    rownames(mat_coloc) <- colnames(m)
    colnames(mat_coloc) <- colnames(m)
    rownames(mat_avoid) <- colnames(m)
    colnames(mat_avoid) <- colnames(m)
    mat_loc <- matrix(ncol = ncol(m), nrow = ncol(m))
    rownames(mat_loc) <- colnames(m)
    colnames(mat_loc) <- colnames(m)

    presence <- presence(spe, method)

    # Calculate pairwise p-values
    for (i in 1:length(rownames(mat_coloc))) {
      for (j in 1:length(colnames(mat_coloc))) {
        p_values <- cell_pair_localization(spe, method = method, cell_type_1 = rownames(mat_coloc)[i], cell_type_2 = rownames(mat_coloc)[j], density = F, distance = distance, presence = presence)
        p_value_coloc <- p_values["coloc_p"]
        mat_coloc[i, j] <- p_value_coloc
        p_value_avoid <- p_values["avoid_p"]
        mat_avoid[i, j] <- p_value_avoid
      }

        # creat coloc and avoid matrix
        mat_loc[upper.tri(mat_loc)] <- mat_avoid[upper.tri(mat_avoid)]
        mat_loc[lower.tri(mat_loc)] <- mat_coloc[lower.tri(mat_coloc)]

        # Select colours
        mycols <- circlize::colorRamp2(
          breaks = c(0, alpha, 1),
          colors = c("red", "white", "black"))

        # Create Heatmap
       ht1 <-  Heatmap(mat_loc,
                heatmap_legend_param = list(
                  title = "Colocalization", at = c(0, alpha-2*alpha/3, alpha-alpha/3 ,alpha, 1),
                  break_dist = 1,
                  labels = c("***", "**", "*", alpha, "ns")),
                  cluster_rows = F, cluster_columns = F, col = mycols
        )
       draw(ht1)
    }
  }

## use corrplot and p-value representation
  if (correlation == TRUE) {
    correlation <- corr.test(m, adjust = "none")$r
    cor_prob <- corr.test(m, adjust = "none")$p

    ht3 <-Heatmap(correlation,
                  rect_gp = gpar(type = "none"), cluster_rows = F, cluster_columns = F, name = "Correlation",
                  cell_fun = function(j, i, x, y, w, h, fill) {
                    if (i > j) {
                      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                    }
                    else if(i<j){
                      grid.text(round(correlation[i,j], digits = 3), x ,y, gp = gpar(fontsize = 10))
                    }
                  }
                 )
    ht4 <-Heatmap(cor_prob,
                  heatmap_legend_param = list(
                    title = "Correlation probability", at = c(0, alpha-2*alpha/3, alpha-alpha/3 ,alpha, 1),
                    break_dist = 1,
                    labels = c("***", "**", "*", alpha, "ns")),
                  rect_gp = gpar(type = "none"), cluster_rows = F, cluster_columns = F, name = "Correlation probability",
                  cell_fun = function(j, i, x, y, w, h, fill) {
                    if (i >j) {
                      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                    }
                    else if(i<j){
                      grid.text(round(cor_prob[i,j], digits = 3), x ,y, gp = gpar(fontsize = 10))
                    }
                  }
    )

    draw(ht3, ht_gap = unit(-75, "mm"))
    draw(ht4, ht_gap = unit(-75, "mm"))

    if (matrix == TRUE) {
      print("correlation")
      print(correlation)
      print("Correlation probability")
      print(cor_prob)
    }
  }
}
