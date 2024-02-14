#' Heatmap correlation of cell types based on gene expression values
#'
#' Generates a heatmap of cell type correlations and p-values as asterisks based on signature matrix
#'
#' @param sig signature matrix
#' @param cor_method correlation method "pearson" or "spearman"
#' @param log logarithmic transformation of signature matrix
#' @param matrix return correlation and correlation probability values as matrix
#'
#' @returns correlation heatmap of cell types based on gene expression
#'
#' @export

corr_expr <- function(sig, log = FALSE, cor_method = "pearson", matrix = FALSE) {
  if (!cor_method %in% c("pearson", "spearman")) {
    stop("cor_method must be either 'pearson' or 'spearman'")
  }

  # option to log scale
  if (log == TRUE) {
    if (any(sig <= 0)) {
      warning("Negative or zero values detected, which are incompatible with log transformation. Adding 1 to all elements.")
    }
    sig <- log(sig + 1)
  }

  # get correlation and p value
  corr_results <- corr.test(sig, method = cor_method, adjust = "none")
  corr_matrix <- corr_results$r
  corr_p <- corr_results$p

  # plot
  corrplot(corr_matrix,
    p.mat = corr_p, method = "color", diag = FALSE, type = "lower",
    sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
    insig = "label_sig", pch.col = "black", tl.col = "black",
    title = "Correlation expression with significance",
    mar = c(0, 0, 2, 0),
    cl.cex = 1.2,
    tl.cex = 1.2,
    cex.main = 1.3,
  )
  cat(paste0("\t", c("Significance correlation:", "* <0.05", "** <0.01", "*** <0.001"), "\n"))

  # optional matrix output
  if (matrix == TRUE) {
    print("Correlation of expression")
    print(corr_matrix)
    print("p_correlation")
    print(corr_p)
  }
}
