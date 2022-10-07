#' Plot Spatial Object
#'
#' @param spatial_obj SpatialExperiment containing deconvolution results
#' @param sample sample to plot
#' @param cell_type Cell Type to plot
#' @param plot_density (default = FALSE) wheter to plot the density
#' @param spot_size size of the spots in the plot
#' @param show_image whether to show the histology image in the background
#' @export
plot_celltype <- function(spatial_obj, sample = "sample01", cell_type = NULL, plot_density = TRUE, spot_size = 1, show_image = TRUE) {
  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is missing or null, but is required")
  }

  if (is.null(cell_type)) {
    stop("Parameter 'cell_types' is null or missing, but is required")
  }

  # check if sample exists in colData
  if (!(sample %in% unlist(SingleCellExperiment::colData(spatial_obj)$sample_id))) {
    stop("Provided sample name not present in object")
  }

  # check if cell types exist in colData
  if (!(cell_type %in% names(SingleCellExperiment::colData(spatial_obj)))) {
    stop("Provided cell types are not available in spatial object")
  }

  # make plots
  spatial <- spatialLIBD::vis_gene(
    spe = spatial_obj,
    sampleid = sample,
    geneid = cell_type,
    point_size = spot_size,
    spatial = show_image
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
#' @param show_image whether to show the histology image in the background
#' @export
plot_cells_per_spot <- function(spatial_obj, plot_type = "spatial", threshold = 0, spot_size = 1.5, show_image = TRUE) {
  if (is.null(spatial_obj)) {
    stop("Paramter 'spatial_obj' is missing or null, but is required")
  }

  if (!plot_type %in% c("spatial", "bar")) {
    stop("Plot_type not supported")
  }

  mat <- get_results_from_object(spatial_obj)


  # initialize result matrix with zeros
  res <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))

  # single threshold or vector
  if (length(threshold) == 1) {
    res[mat > threshold] <- 1
  } else if (length(threshold) == nrow(mat)) {
    # vector matches matrix
    for (i in 1:nrow(mat)) {
      res[i, mat[i, ] > threshold[i]] <- 1
    }
  }

  rownames(res) <- rownames(mat)
  colnames(res) <- colnames(mat)

  # count cells
  plot_data <- data.frame(spot = rownames(res), value = as.integer(apply(res, 1, sum)))

  # make plot
  if (plot_type == "bar") {
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
  } else if (plot_type == "spatial") {
    obj <- spatial_obj
    tmp <- SingleCellExperiment::colData(obj)
    SummarizedExperiment::colData(obj) <- cbind(tmp, value = plot_data[["value"]])

    pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))

    plot <- spatialLIBD::vis_clus(
      obj,
      sampleid = "sample01",
      clustervar = "value",
      point_size = spot_size,
      spatial = show_image,
      colors = rev(pal(ncol(tmp)))
    )
  } else {
    plot <- NULL
  }
  plot
}

#' Function to plot deconvolution results
#'
#' Generate Hex Plot of a SpatialExperiment containing deconvolution results
#'
#' @param spe deconvolution result in Form of a SpatialExperiment
#' @param cell_type one or more celltype to plot
#' @param sample_id sample id to plot, default: "sample01"
#' @param image_id which image to plot, default: "lowres"
#'
#' @returns plot of cell type fractions
#'
#' @export
#' @example
#' # TODO
new_plot_celltype <- function(spe, cell_type = NULL, sample_id = "sample01", image_id = "lowres") {
  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if (is.null(cell_type)) {
    stop("Parameter 'cell_type' is null or missing, but is required")
  }

  # check that celltypes are present in object
  if (!all(cell_type %in% names(colData(spe)))) {
    stop("Provides cell types are not presend in SpatialExperiment")
  }

  df <- as.data.frame(cbind(SpatialExperiment::spatialCoords(spe), colData(spe)))

  # img = SpatialExperiment::imgRaster(spe), for later

  # scale coordinates with scalefactor
  df$pxl_col_in_fullres <- df$pxl_col_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sample_id, image_id = image_id)
  df$pxl_row_in_fullres <- df$pxl_row_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sample_id, image_id = image_id)

  # need to work on this, does this work when first spot is completely separated from the rest?
  distance_guess <- min(sqrt((df$pxl_col_in_fullres[1] - df$pxl_col_in_fullres[-1])^2 + (df$pxl_row_in_fullres[1] - df$pxl_row_in_fullres[-1])^2))

  # build_hex coordinates
  get_hex_polygon <- function(x, y, dist) {
    angle <- seq(0, 2 * pi, length.out = 7)[-7] # angles of
    res <- cbind(
      x = x + sin(angle + 0.5) * dist / 2, # rotated by 30 degrees
      y = y + cos(angle + 0.5) * dist / 2
    )
    return(rbind(res, res[1, ]))
  }

  # get coordinates for all hex polygons
  get_polygon_geometry <- function(grid, dist) {
    res <- list()
    for (i in seq_len(nrow(grid))) {
      res <- c(res, list(sf::st_polygon(list(get_hex_polygon(grid$pxl_col_in_fullres[i], grid$pxl_row_in_fullres[i], dist)))))
    }
    return(sf::st_sfc(res)) # Convert to 'simple feature collection'
  }

  # polygons
  new_geom <- get_polygon_geometry(df, distance_guess)

  # preparing the dataframe with sf
  sf_points <- sf::st_as_sf(df, coords = c("pxl_col_in_fullres", "pxl_row_in_fullres"))

  # no overwrite the polygon points with hex polygons
  sf_poly <- sf::st_set_geometry(sf_points, new_geom)

  # plot(sf_poly[cell_type], border = NA)
  p <- ggplot() +
    geom_sf(aes_string(fill = cell_type), lwd = 0, data = sf_poly) +
    colorspace::scale_fill_continuous_sequential("Rocket") +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(), panel.grid = element_blank()
    )

  return(p)

  # TODO
  # add facet wrap
  # add color selection
  # add background image
  # try to reducet the gap
}
