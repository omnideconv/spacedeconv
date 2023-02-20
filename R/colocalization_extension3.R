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
        mat_loc[upper.tri(mat_loc, diag = T)] <- mat_avoid[upper.tri(mat_avoid, diag = T)]
        mat_loc[lower.tri(mat_loc)] <- mat_coloc[lower.tri(mat_coloc)]


        # Select colours
        mycols <- circlize::colorRamp2(
          breaks = c(0, 0.05, 1),
          colors = c("red","white", "green"))


         # Create Heatmap


        cell_fun = function(j, i, x, y, w, h, fill){
        if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
        grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }
        if (mat_loc[i, j]  < 0.01 & as.numeric(x) <= 1 - as.numeric(y) + 1e-6){
        grid.text(paste0(sprintf("%.2f", mat_loc[i, j]),"**"), x, y, gp = gpar(fontsize = 10))
        } else if (mat_loc[i, j]  <= 0.05 & as.numeric(x) <= 1 - as.numeric(y) + 1e-6){
grid.text(paste0(sprintf("%.2f", mat_loc[i, j]),"*"), x, y, gp = gpar(fontsize = 10))
        }
          grid.text(paste0(sprintf("%.2f", mat_loc[i, j]),"*"), x, y, gp = gpar(fontsize = 10))
        }
      }

       ht1 <-  Heatmap(mat_loc,

                heatmap_legend_param = list(
                  title = "Colocalization/Avoidance", at = c(0.001, 0.01, 0.05, 1),
                  break_dist = 1,
                  labels = c("0.001", "0.01", "0.05", "ns")),
                  cluster_rows = F, cluster_columns = F, col = mycols,
                  column_title = "Avoidance", row_title = "Colocalization",
                  cell_fun = cell_fun,

        )
       draw(ht1)
    }
  }

  if (correlation == TRUE) {
    correlation <- corr.test(m, adjust = "none")$r
    cor_prob <- corr.test(m, adjust = "none")$p
    corrplot(correlation, p.mat = cor_prob, method = 'color', diag = FALSE, type = 'lower',
             sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
             insig = 'label_sig', pch.col = 'black', tl.col = "black",
            title = "Correlation with significance",
            mar=c(0,0,2,0))
    cat(paste0("\t", c("Significance:", "* <0.05", "** <0.01", "*** <0.001"), "\n"))

  }
}
