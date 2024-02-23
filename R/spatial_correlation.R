#library(spacedeconv)
# WARNING: the following functions were not found although spacedeconv was loaded
# "cor.mtest"

#library(corrplot)

#' create a correlation plot
#' @param spe the spe object
#' @param method method to used for correlation analysis - deconvolution method or decoupleR e.g. "cell2location"
#' @param adjust method used to adjust p-values for multiple testing, see stats::p.asjust.methods for details
#' @param variables if not provided, the function selects variables based on the specified method, it can also be a vector
#' @param sig.level it can also be a vector
#' @param type Character, 'full' 'upper' or 'lower' (default), display full matrix, lower triangular or upper triangular matrix.
#' @param diag Logical, whether display the correlation coefficients on the principal diagonal. default is FALSE
#' @param order Character, the ordering method of the correlation matrix.'original' for original order (default),'AOE' for the angular order of the eigenvectors, 'FPC' for the first principal component order, 'hclust' for the hierarchical clustering order,'alphabet' for alphabetical order.
#' @param insig If 'blank', wipe away the corresponding glyphs; if 'p-value', add p-values the corresponding glyphs; if 'pch', add characters (see pch for details) on corresponding glyphs; if 'n', don't take any measures; if 'label_sig' (default), mark significant correlations with pch (see sig.level).
#' @param plot_layout represents the method parameter in the original corrplot function.The layout of the correlation plot. color is default. Choose between "circle", "square", "ellipse", "number", "shade", "color", "pie"
#' @param addCoef.col to add labels to the plot showing the correlation values, default is NULL, the user can choose any other color to add and color the values.
#' @param ... additional parameters passed to corrplot function
#' @export

spatialcorr <- function(spe,
                        method, # method to used for correlation analysis - deconvolution method or decoupleR e.g. "cell2location"
                        adjust = "fdr", # method used to adjust p-values for multiple testing
                        variables = NULL, # if not provided, the function selects variables based on the specified method
                        sig.level = 0.05,
                        type = "lower",
                        diag = FALSE,
                        order = "original",
                        plot_layout = "color",
                        insig = "label_sig", # to label or not the significant values
                        addCoef.col = NULL,
                        ...) {
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

  # check if valid type
  valid_type <- c("full", "lower", "upper")
  if (!(type %in% valid_type)) {
    stop("Invalid argument provided for the 'type' parameter. Choose one of: ", paste(valid_type, collapse = ", "))
  }

  # Check if diag is logical
  if (!is.logical(diag)) {
    stop("Parameter 'diag' must be logical (TRUE or FALSE).")
  }

  # check if valid order
  valid_order <- c("original", "AOE", "FPC", "hclust", "alphabet")
  if (!(order %in% valid_order)) {
    stop("Invalid argument provided for the 'order' parameter. Choose one of: ", paste(valid_order, collapse = ", "))
  }

  # check if valid plot_layout
  valid_plot <- c("circle", "square", "ellipse", "number", "shade", "color", "pie")
  if (!(plot_layout %in% valid_plot)) {
    stop("Invalid argument provided for the 'order' parameter. Choose one of: ", paste(valid_plot, collapse = ", "))
  }

  # check if valid plot_layout
  if (plot_layout == "number") {
    print("Consider setting the 'insig' parameter to: 'blank', otherwise the significance markers (asterisks or X) would overlap with the correlation values")
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
    method = plot_layout, # default is color
    diag = diag, # default is FALSE, or defined above
    type = type, # default is "lower", or else defined above
    sig.level = sig.level,
    insig = insig, # add significant level asterisks
    pch.cex = 0.9,
    pch.col = "black",
    tl.col = "black",
    mar = c(0, 0, 2, 0),
    cex.lab = 1.6,
    cex.main = 1.5,
    order = order, # default is "original", or else defined above
    addCoef.col = addCoef.col, # default is NULL
    ...
  )

  # Return matrix of correlations and adjusted p-values
  invisible(list(corr = r, padj = p.adj))
}
