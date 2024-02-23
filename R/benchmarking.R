#' Benchmarking scatterplot to compare two spatial objects
#'
#' @param spe1 SpatialExperiment
#' @param value1 deconvolution result to plot
#' @param spe2 SpatialExperiment
#' @param value2 deconvolution result to plot
#'
#' @export
plot_scatter <- function(spe1, value1, spe2, value2, log_scale = FALSE) {
  # check object class
  if (!is(spe1, "SpatialExperiment") || !is(spe2, "SpatialExperiment")) {
    cli::cli_alert_danger("Provided objects have to be SpatialExperiments")
    stop()
  }

  # check value availability
  if (!checkCol(spe1, value1)) {
    cli::cli_alert_danger("Provided Column for object 1 is not available")
    stop()
  }

  if (!checkCol(spe2, value2)) {
    cli::cli_alert_danger("Provided Column for object 2 is not available")
    stop()
  }

  # check that both objects have the same number of spots
  if (ncol(spe1) != ncol(spe2)) {
    cli::cli_alert_warning("Spatial Objects have different number of spots")
  }

  # construct data frame for plotting
  df1 <- colData(spe1)[, value1, drop = FALSE]
  colnames(df1) <- "value1"
  df1$spot <- rownames(df1)

  df2 <- colData(spe2)[, value2, drop = FALSE]
  colnames(df2) <- "value2"
  df2$spot <- rownames(df2)

  df <- merge(df1, df2, by = "spot")

  cor_value <- cor(df$value1, df$value2, use = "complete.obs") # handle NA values

  # construct plot
  plot <- ggplot(df, aes(x = value1, y = value2)) +
    geom_point(shape = 21, color = "blue", fill = "lightblue", size = 1, stroke = 0.5) +
    geom_abline(slope = 1, linetype = "dashed", col="red") +
    xlab(paste(value1)) +
    ylab(paste(value2)) +
    coord_fixed(ratio = 1) +
    geom_text(
      x = Inf, y = Inf, label = paste("Correlation:", round(cor_value, 2)),
      hjust = 1.1, vjust = 1.1, size = 5
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    ) +
    ggtitle(paste(value1, "vs.", value2))


  # Apply log scale if log_scale is TRUE
  if (log_scale) {
    plot <- plot + scale_x_log10() + scale_y_log10()
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
