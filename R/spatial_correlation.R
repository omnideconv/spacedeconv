# library(spacedeconv)
# WARNING: the following functions were not found although spacedeconv was loaded
# "cor.mtest"

# library(corrplot)

#' create a correlation plot
#' @param spe #the spe object
#' @param method #method to used for correlation analysis - deconvolution method or decoupleR e.g. "cell2location"
#' @param adjust #method used to adjust p-values for multiple testing
#' @param variables #if not provided, the function selects variables based on the specified method, it can also be a vector
#' @param sig.level #it can also be a vector
#' @export

spatialcorr <- function(spe,
                        method, # method to used for correlation analysis - deconvolution method or decoupleR e.g. "cell2location"
                        adjust = "fdr", # method used to adjust p-values for multiple testing
                        variables = NULL, # if not provided, the function selects variables based on the specified method
                        sig.level = 0.05) {
  # CHECKS
  # spe
  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required.")
  }

  # method availability in the spe object
  if (!any(startsWith(method, available_results(spe)) | startsWith(available_results(spe), method))) {
    stop("The specified 'method' does not exist in the spe.")
  }

  # method and variables missing
  if (is.null(method) && is.null(variables)) {
    stop("Parameter 'method' and 'variables' is null or missing, but one of them is required.")
  }

  # adjust - check if its a valid method
  valid_adjust_methods <- c("bonferroni", "holm", "hochberg", "hommel", "fdr")
  if (!(adjust %in% valid_adjust_methods)) {
    stop("Invalid adjustment method. Choose one of: ", paste(valid_adjust_methods, collapse = ", "))
  }

  # variables - check if its a vector and if they exist in the object
  if (!is.null(variables) && !is.vector(variables)) {
    stop("'variables' must be a vector.")
  }

  if (!is.null(variables) && any(!variables %in% colnames(colData(spe)))) {
    stop("One or more specified variables do not exist in the spe object.")
  }

  # sig.level
  if (sig.level <= 0 || sig.level >= 1) {
    stop("'sig.level' must be between 0 and 1.")
  }

  # Select variables of interests
  # If variables is not provided but method is specified, the function selects variables based on those available in the dataset that start with the specified method.
  if (is.null(variables) && !is.null(method)) {
    selvar <- available_results(spe)[startsWith(available_results(spe), method)]

    # If variables is provided, it selects the intersection of the specified variables and the available variables in the dataset.
  } else {
    selvar <- intersect(available_results(spe), variables)
  }
  m <- as.matrix(colData(spe)[, selvar])

  cor_res <- cor.mtest(m)
  r <- cor_res$uppCI # correlation
  p <- cor_res$p # p-value

  # Correct for multiple testing
  tmp <- p.adjust(p[lower.tri(p)], method = adjust) # Adjust p-values for multiple testing using the specified method
  p.adj <- matrix(1, ncol = ncol(p), nrow = nrow(p))
  colnames(p.adj) <- colnames(r)
  rownames(p.adj) <- rownames(r)
  p.adj[lower.tri(p.adj)] <- tmp
  p.adj[upper.tri(p.adj)] <- t(p.adj)[upper.tri(p.adj)] # symmetrical matrix
  diag(p.adj) <- 0 # diag set to 0

  corrplot(r,
    p.mat = p.adj,
    method = "color",
    diag = FALSE,
    type = "lower", # the plot contains only the lower part
    sig.level = sig.level,
    insig = "label_sig", # add significant level asterisks
    pch.cex = 0.9,
    pch.col = "black",
    tl.col = "black",
    mar = c(0, 0, 2, 0),
    cex.lab = 1.6,
    cex.main = 1.5
  )

  # Return matrix of correlations and adjusted p-values
  invisible(list(corr = r, padj = p.adj))
}
