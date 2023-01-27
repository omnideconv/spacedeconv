#' Heatmaps cell pair localization and correlation of scores
#'
#' Generates a heatmap of cell pair localization p-values and a correlation heatmap of scores between all celltypes found after deconvolution
#'
#' @param spe SpatialExperiment object
#' @param method deconvolution method
#' @param distance size area of interest
#' @param correlation correlation heatmap based on scores
#' @param localization colocalization and avoidance heatmap
#'
#' @returns cell pair localization heatmap and correlation heatmap
#'
#' @export

localization_heatmap <- function(spe, method, distance = 0, correlation = TRUE, localization = TRUE, matrix = FALSE) {
  # create matrix with scores for each spot and celltype
  available <- available_results(spe)[startsWith(available_results(spe), method)]
  m <- as.matrix(colData(spe)[, available])

  if (localization == TRUE) {
    # Create empty matrix to fill in p-values
    mat_coloc <- matrix(ncol = ncol(m), nrow = ncol(m))
    rownames(mat_coloc) <- colnames(m)
    colnames(mat_coloc) <- colnames(m)
    mat_avoid <- matrix(ncol = ncol(m), nrow = ncol(m))
    rownames(mat_avoid) <- colnames(m)
    colnames(mat_avoid) <- colnames(m)

    presence <- presence(spe, method)

    # Calculate pairwise p-values
    for (i in 1:length(rownames(mat_coloc))) {
      for (j in 1:length(colnames(mat_coloc))) {
        p_values <- cell_pair_localization(spe, method = method, cell_type_1 = rownames(mat_coloc)[i], cell_type_2 = rownames(mat_coloc)[j], density = F, distance = distance, presence = presence)
        p_value_coloc <- p_values["coloc_p"]
        mat_coloc[i, j] <- p_value_coloc
        p_value_avoid <- p_values["avoid_p"]
        mat_avoid[i, j] <- p_value_avoid

        # Select colours
        mycols <- circlize::colorRamp2(
          breaks = c(0, 0.05, 1),
          colors = c("red", "white", "black")
        )
        # Create heatmap
        ht1 <- Heatmap(mat_coloc,
          rect_gp = gpar(type = "none"), cluster_rows = F, cluster_columns = F, name = "Colocalization significance (left)", col = mycols,
          cell_fun = function(j, i, x, y, w, h, fill) {
            if (i > j) {
              grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
            }
          }
        )

        colnames(mat_avoid) <- NULL
        ht2 <- Heatmap(mat_avoid,
          rect_gp = gpar(type = "none"), cluster_rows = F, cluster_columns = F, name = "Avoidance significance (right)", col = mycols,
          cell_fun = function(j, i, x, y, w, h, fill) {
            if (i < j) {
              grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
            }
          }
        )
      }
    }
    draw(ht1 + ht2, ht_gap = unit(-75, "mm"))
    if (matrix == TRUE) {
      print("mat_coloc")
      print(mat_coloc)
      print("mat_avoid")
      print(mat_avoid)
    }
  }
  if (correlation == TRUE) {
    correlation <- corr.test(m, adjust = "none")$r
    cor_prob <- corr.test(m, adjust = "none")$p
    colnames(cor_prob) <- NULL

    ht3 <-Heatmap(correlation,
                  rect_gp = gpar(type = "none"), cluster_rows = F, cluster_columns = F, name = "Correlation (right)",
                  cell_fun = function(j, i, x, y, w, h, fill) {
                    if (i > j) {
                      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                    }
                  }
                )
    ht4 <-Heatmap(cor_prob,
                  rect_gp = gpar(type = "none"), cluster_rows = F, cluster_columns = F, name = "Correlation probability (left)",
                  cell_fun = function(j, i, x, y, w, h, fill) {
                    if (i < j) {
                      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                    }
                  }
                 )

    draw(ht3 + ht4, ht_gap = unit(-75, "mm"))

    if (matrix == TRUE) {
      print("correlation")
      print(correlation)
      print("Correlation probability")
      print(cor_probability)
    }
  }
}
