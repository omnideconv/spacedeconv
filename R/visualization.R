#' Plot Spatial Object
#'
#' @param spatial_obj SpatialExperiment containing deconvolution results
#' @param sample sample to plot
#' @param cell_type Cell Type to plot
#' @export
plot_celltype <- function(spatial_obj, sample = "sample01", cell_type = NULL) {
  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is missing or null, but is required")
  }

  if (is.null(cell_type)) {
    stop("Parameter 'cell_types' is null or missing, but is required")
  }

  # check if cell types exist in colData
  if (!(cell_type %in% names(SingleCellExperiment::colData(spatial_obj)))) {
    stop("provided cell types are not available in spatial object")
  }

  # make plots
  spatial <- spatialLIBD::vis_gene(
    spe = spatial_obj, sampleid = sample,
    geneid = cell_type, point_size = 2
  )

  # extract distribution from object

  data <- data.frame(values = SingleCellExperiment::colData(spatial_obj)[[cell_type]])

  density <- ggplot2::ggplot(data, mapping = ggplot2::aes_string(x = "values", fill = "cell_type")) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = mean(unlist(data))),
      color = "blue",
      linetype = "dashed",
      size = 1
    ) +
    ggplot2::theme_classic()

  cowplot::plot_grid(spatial, density, labels = "AUTO", scale = c(1, 0.85))
}
