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

corr_expr <- function(sig, log = FALSE, cor_method = c("pearson", "spearman"), matrix = FALSE) {
  if (log == TRUE) {
    sig <- log(sig + 1)
  }
  corr_matrix <- corr.test(sig, method = cor_method, adjust = "none")$r
  corr_p <- corr.test(sig, adjust = "none")$p
  corrplot(corr_matrix,
    p.mat = corr_p, method = "color", diag = FALSE, type = "lower",
    sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
    insig = "label_sig", pch.col = "black", tl.col = "black",
    title = "Correlation expression with significance",
    mar = c(0, 0, 2, 0)
  )
  cat(paste0("\t", c("Significance correlation:", "* <0.05", "** <0.01", "*** <0.001"), "\n"))
  if (matrix == TRUE) {
    print("Correlation of expression")
    print(corr_matrix)
    print("p_correlation")
    print(corr_p)
  }
}
