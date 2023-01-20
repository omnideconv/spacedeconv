#' Benchmarking scatterplot to compare two spatial objects
#'
#' @param spe1 SpatialExperiment
#' @param value1 deconvolution result to plot
#' @param spe2 SpatialExperiment
#' @param value2 deconvolution result to plot
#'
#' @export
plot_scatter <- function(spe1, value1, spe2, value2) {
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

  # construct data frame for plotting
  df1 <- colData(spe1)[value1]
  colnames(df1) <- "value1"
  df1$spot <- rownames(df1)

  df2 <- colData(spe2)[value2]
  colnames(df2) <- "value2"
  df2$spot <- rownames(df2)

  df <- merge(df1, df2, by = "spot")


  # construct plot
  plot <- ggplot(data.frame(df), aes(x = value1, y = value2)) +
    geom_point() +
    geom_abline(slope = 1) +
    xlab(value1) +
    ylab(value2)

  return(plot)
}
