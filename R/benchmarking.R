#' Compare Two Result Columns with a Scatterplot
#'
#' Creates a scatterplot to compare two result columns, either from a single
#' `SpatialExperiment` (use `spe` + `value1` + `value2`) or from two different
#' objects (use `spe1` + `spe2` + `value1` + `value2`).
#'
#' @param spe SpatialExperiment with quantification results
#' @param value1 Numeric colData column to plot on the x-axis (e.g. cell abundance or activity).
#' @param value2 Numeric colData column to plot on the y-axis.
#' @param spe1 First SpatialExperiment object.
#' @param spe2 Second SpatialExperiment object.
#' @param log_scale Logical, whether to use log scale for axes.
#' @param dot_color Color for the dots in the plot.
#' @param point_alpha Alpha for the dots, controlling transparency.
#' @param fix_coords Logical, whether to use equal axis ranges and a 1:1 aspect ratio.
#' @param coord_range Optional increasing numeric vector of length 2 setting both
#'   axis limits in original data units. NULL uses automatic limits.
#' @param title Title of the plot, or NULL to omit it.
#' @param xlab,ylab Axis labels; default to the selected column names. Use NULL
#'   to omit a label.
#'
#' @param point_size Size of the points, in ggplot2 size units.
#' @param stats_size Font size of the Pearson r and p-value annotation, in points.
#'
#' @details The annotation reports Pearson's r and its standard correlation-test
#'   p-value, which does not account for spatial dependence between spots.
#'   The solid line is a linear regression; the dashed line represents equality.
#'
#' @return A scatterplot.
#' @export
plot_scatter <- function(spe = NULL, value1, value2, spe1 = NULL, spe2 = NULL, log_scale = FALSE, dot_color = "#1f77b4", point_alpha = 0.8, fix_coords = FALSE, coord_range = NULL, title = "Comparative Scatterplot", xlab = value1, ylab = value2, point_size = 2, stats_size = 15) {
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

  if (!is.null(coord_range) &&
      (!is.numeric(coord_range) || length(coord_range) != 2L ||
       any(!is.finite(coord_range)) || coord_range[1] >= coord_range[2])) {
    stop("coord_range must contain two finite, increasing numeric limits.")
  }
  if (log_scale && !is.null(coord_range) && any(coord_range <= 0)) {
    stop("coord_range must be positive when log_scale = TRUE.")
  }

  # Construct plot
  plot <- ggplot(df, aes(x = value1, y = value2)) +
    geom_point(color = dot_color, size = point_size, alpha = point_alpha, shape = 16) +
    ggplot2::geom_smooth(method = "lm", color = "#454545", se = FALSE) + # Regression line
    geom_abline(slope = 1, linetype = "dashed", col = "gray") +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 21, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(size = 17),
      axis.text = ggplot2::element_text(size = 14),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey92")
    ) +
    ggplot2::ggtitle(title)

  p_label <- format.pval(p_value, digits = 2, eps = 2e-16)
  p_label <- if (startsWith(p_label, "<")) paste("p", p_label) else paste("p =", p_label)
  # Anchor inside the panel, independently of data limits and axis transforms.
  plot <- plot + ggplot2::annotation_custom(grid::textGrob(
    label = paste0(sprintf("r = %.2f", cor_value), "\n", p_label),
    x = grid::unit(0.05, "npc"), y = grid::unit(0.95, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(col = "#454545", fontsize = stats_size)))

  # Apply log scale if log_scale is TRUE
  if (log_scale) {
    plot <- plot + ggplot2::scale_x_log10() + ggplot2::scale_y_log10()
  }

  # Share ranges when requested, retaining negative values for activity scores.
  limits <- coord_range
  if (fix_coords && is.null(limits)) {
    values <- c(df$value1, df$value2)
    values <- values[is.finite(values)]
    if (log_scale) values <- values[values > 0]
    limits <- range(values)
  }
  if (fix_coords) {
    plot <- plot + ggplot2::coord_fixed(ratio = 1, xlim = limits, ylim = limits)
  } else if (!is.null(limits)) {
    plot <- plot + ggplot2::coord_cartesian(xlim = limits, ylim = limits)
  }

  return(plot)
}





#' Compare Signatures
#'
#' Creates a scatterplot comparing two signature matrices across shared genes.
#'
#' @param signature1 Signature matrix.
#' @param signature2 Signature matrix.
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
