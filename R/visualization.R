#' Plot Spatial Object
#'
#' @param spatial_obj SpatialExperiment containing deconvolution results
#' @param sample sample to plot
#' @param cell_type Cell Type to plot
#' @param plot_density (default = FALSE) wheter to plot the density
#' @param spot_size size of the spots in the plot
#' @export
plot_celltype <- function(spatial_obj, sample = "sample01", cell_type = NULL, plot_density = TRUE, spot_size = 1) {
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
    geneid = cell_type, point_size = spot_size
  )

  plot <- spatial

  if (plot_density) {
    # extract distribution from object
    data <- data.frame(values = SingleCellExperiment::colData(spatial_obj)[[cell_type]], id = rep(cell_type, nrow(SingleCellExperiment::colData(spatial_obj))))
    density <- ggplot2::ggplot(data, mapping = ggplot2::aes_string(x = "values", y = "id")) + # fill... see ggridges docs
      # ggplot2::geom_density() +
      ggridges::geom_density_ridges() + # _gradient()
      # ggplot2::scale_fill_viridis_c() +
      ggplot2::scale_y_discrete() +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = mean(unlist(data["values"]))),
        color = "red",
        linetype = "dashed",
        size = 1
      ) +
      ggplot2::theme_classic() +
      # ggplot2::ylim(c(0, 1000)) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank()
      )
    # ggplot2::ylim(0, max(data["values"]))

    # cowplot::plot_grid(spatial, density)
    plot <- gridExtra::grid.arrange(spatial, density, ncol = 2)
  }

  plot
}

#' Plot Cells per Spot
#'
#' @param spatial_obj SpatialExperiment containing deconvolution results
#' @param plot_type bar chart or spatial (bar, spatial)
#' @param threshold threshold for presence/absence, single value or vector of length nrow(spatial_obj)
#' @param spot_size size of the dots
#' @export
plot_cells_per_spot <- function(spatial_obj, plot_type = "spatial", threshold = 0, spot_size = 1.5) {
  if (is.null(spatial_obj)) {
    stop("Paramter 'spatial_obj' is missing or null, but is required")
  }

  mat <- get_results_from_object(spatial_obj)


  # initialize result matrix with zeros
  res <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))

  # single threshold or vector
  if (length(threshold) == 1) {
    res[mat > threshold] <- 1
  } else if (length(threshold) == nrow(mat)) { # vector matches matrix
    for (i in 1:nrow(mat)) {
      res[i, mat[i, ] > threshold[i]] <- 1
    }
  }

  rownames(res) <- rownames(mat)
  colnames(res) <- colnames(mat)

  # count cells
  plot_data <- data.frame(spot = rownames(res), value = as.integer(apply(res, 1, sum)))

  # make plot
  if (plot_type=="bar"){
    plot <- ggplot2::ggplot(plot_data) +
      ggplot2::geom_bar(aes_string(x = "value"), stat = "count", fill = "black") +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = mean(unlist(plot_data["value"]))),
        color = "red",
        linetype = "dashed",
        size = 1
      ) +
      ggplot2::geom_hline(yintercept = 0, size = 1.2) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank(),
      )
  } else if (plot_type=="spatial"){
    obj = spatial_obj
    tmp = SingleCellExperiment::colData(obj)
    SingleCellExperiment::colData(obj) <- cbind(tmp, value=plot_data[["value"]])

    pal = colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))

    plot <- spatialLIBD::vis_clus(obj, sampleid = "sample01", clustervar ="value", point_size = spot_size, colors = rev(pal(ncol(tmp))))
  } else {
    plot = NULL
  }
  plot
}
