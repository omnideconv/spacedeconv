#' Comparative Scatterplot of Two SpatialExperiments
#'
#' Creates a scatterplot to compare deconvolution results between two SpatialExperiment objects.
#'
#' @param spe SpatialExperiment with quantification results
#' @param value1 Deconvolution result to plot from the first object.
#' @param value2 Deconvolution result to plot from the second object.
#' @param spe1 First SpatialExperiment object.
#' @param spe2 Second SpatialExperiment object.
#' @param log_scale Logical, whether to use log scale for axes.
#' @param dot_color Color for the dots in the plot.
#' @param point_alpha Alpha for the dots, controlling transparency.
#' @param fix_coords Logical, whether to fix coordinate system to be equal.
#' @param coord_range Numeric vector of length 2 to set coordinate limits.
#' @param title Title of the plot.
#'
#' @return A scatterplot
#' @export
plot_scatter <- function(spe = NULL, value1, value2, spe1 = NULL, spe2 = NULL, log_scale = FALSE, dot_color = "#1f77b4", point_alpha = 0.8, fix_coords = FALSE, coord_range = NULL, title = "Comparative Scatterplot") {
  if (!is.null(spe)) {
    if (!is(spe, "SpatialExperiment")) {
      cli::cli_alert_danger("Provided object has to be a SpatialExperiment")
      stop()
    }

    # Check value availability
    if (!checkCol(spe, value1)) {
      cli::cli_alert_danger("Provided Column for value1 is not available in the provided SpatialExperiment")
      stop()
    }

    if (!checkCol(spe, value2)) {
      cli::cli_alert_danger("Provided Column for value2 is not available in the provided SpatialExperiment")
      stop()
    }

    # Construct data frame for plotting
    df <- colData(spe)[, c(value1, value2), drop = FALSE]
    colnames(df) <- c("value1", "value2")
    df$spot <- rownames(df)
  } else {
    if (is.null(spe1) || is.null(spe2)) {
      cli::cli_alert_danger("Both spe1 and spe2 must be provided if spe is not specified")
      stop()
    }

    if (!is(spe1, "SpatialExperiment") || !is(spe2, "SpatialExperiment")) {
      cli::cli_alert_danger("Provided objects have to be SpatialExperiments")
      stop()
    }

    # Check value availability
    if (!checkCol(spe1, value1)) {
      cli::cli_alert_danger("Provided Column for value1 is not available in spe1")
      stop()
    }

    if (!checkCol(spe2, value2)) {
      cli::cli_alert_danger("Provided Column for value2 is not available in spe2")
      stop()
    }

    # Check that both objects have the same number of spots
    if (ncol(spe1) != ncol(spe2)) {
      cli::cli_alert_warning("Spatial Objects have different number of spots")
    }

    # Construct data frame for plotting
    df1 <- colData(spe1)[, value1, drop = FALSE]
    colnames(df1) <- "value1"
    df1$spot <- rownames(df1)

    df2 <- colData(spe2)[, value2, drop = FALSE]
    colnames(df2) <- "value2"
    df2$spot <- rownames(df2)

    df <- merge(df1, df2, by = "spot")
  }

  cor_test_result <- cor.test(df$value1, df$value2, use = "complete.obs") # Handle NA values
  cor_value <- cor_test_result$estimate
  p_value <- cor_test_result$p.value

  # Set coordinate range
  max_value <- max(c(df$value1, df$value2), na.rm = TRUE)
  coord_range <- c(0, max_value + max_value * 0.1)

  # Construct plot
  plot <- ggplot(df, aes(x = value1, y = value2)) +
    geom_point(color = dot_color, size = 2, alpha = point_alpha, shape = 16) +
    geom_smooth(method = "lm", color = "red", se = FALSE) + # Regression line
    geom_abline(slope = 1, linetype = "dashed", col = "gray") +
    xlab(paste(value1)) +
    ylab(paste(value2)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    ) +
    ggtitle(title)

  # Annotate with the correlation coefficient and p-value
  plot <- plot + annotate("text", x = Inf, y = Inf, label = paste0("Corr: ", round(cor_value, 2), ", P val: ", format.pval(p_value, digits = 3)), hjust = 1.1, vjust = 1, color = "red", size = 5)

  # Apply log scale if log_scale is TRUE
  if (log_scale) {
    plot <- plot + scale_x_log10() + scale_y_log10()
  }

  # If fix_coords is TRUE, apply coord_fixed or coord_cartesian based on whether coord_range is specified
  if (fix_coords) {
    if (!is.null(coord_range) && length(coord_range) == 2) {
      plot <- plot + coord_cartesian(xlim = coord_range, ylim = coord_range)
    } else {
      plot <- plot + coord_fixed(ratio = 1)
    }
  }

  return(plot)
}


#' Compare Signatures
#'
#' @param signature1 signature
#' @param signature2 signature
#' @export
compare_signatures <- function(signature1, signature2) {
  df1 <- data.frame(signature1)
  df1$gene <- rownames(signature1)
  df1 <- tidyr::pivot_longer(df1, !gene, values_to = "signature1")

  df2 <- data.frame(signature2)
  df2$gene <- rownames(signature2)
  df2 <- tidyr::pivot_longer(df2, !gene, values_to = "signature2")

  df <- merge(df1, df2)

  plot <- ggplot(df, aes(x = signature1, y = signature2)) +
    geom_point() +
    geom_abline(slope = 1) +
    xlab("signature 1") +
    ylab("signature 2")

  return(plot)
}
