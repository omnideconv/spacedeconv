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

  data <- data.frame(values = SingleCellExperiment::colData(spatial_obj)[[cell_type]], id=rep(cell_type, nrow(SingleCellExperiment::colData(spatial_obj))))
  density <- ggplot2::ggplot(data, mapping = ggplot2::aes_string(x = values, y=id, fill=ggplot2::after_stat(x))) +
    #ggplot2::geom_density() +
    ggridges::geom_density_ridges_gradient() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::scale_y_discrete() +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = mean(unlist(data["values"]))),
      color = "red",
      linetype = "dashed",
      size = 1
    ) +
    ggplot2::theme_classic() +
    #ggplot2::ylim(c(0, 1000)) +
    ggplot2::theme(legend.position = "none",
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank())
    #ggplot2::ylim(0, max(data["values"]))

  cowplot::plot_grid(spatial, density)
}
