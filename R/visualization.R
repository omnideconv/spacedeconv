#   if (plot_density) {
#     # extract distribution from object
#     data <- data.frame(values = SingleCellExperiment::colData(spatial_obj)[[cell_type]], id = rep(cell_type, nrow(SingleCellExperiment::colData(spatial_obj))))
#     density <- ggplot2::ggplot(data, mapping = ggplot2::aes_string(x = "values", y = "id")) + # fill... see ggridges docs
#       # ggplot2::geom_density() +
#       ggridges::geom_density_ridges() + # _gradient()
#       # ggplot2::scale_fill_viridis_c() +
#       ggplot2::scale_y_discrete() +
#       ggplot2::geom_vline(
#         ggplot2::aes(xintercept = mean(unlist(data["values"]))),
#         color = "red",
#         linetype = "dashed",
#         size = 1
#       ) +
#       ggplot2::theme_classic() +
#       # ggplot2::ylim(c(0, 1000)) +
#       ggplot2::theme(
#         legend.position = "none",
#         axis.text.y = ggplot2::element_blank(),
#         axis.ticks.y = ggplot2::element_blank(),
#         axis.line.y = ggplot2::element_blank(),
#         axis.title = ggplot2::element_blank()
#       )
#     # ggplot2::ylim(0, max(data["values"]))
#
#     # cowplot::plot_grid(spatial, density)
#     plot <- gridExtra::grid.arrange(spatial, density, ncol = 2)



#' Plot Cells per Spot
#'
#' @param spatial_obj SpatialExperiment containing deconvolution results
#' @param plot_type bar chart or spatial (bar, spatial)
#' @param threshold threshold for presence/absence, single value or vector of length nrow(spatial_obj)
#' @param sample_id sample of SpatialExperiment to be plotted
#' @param image_id image of SpatialExperiment for the background annotation
#' @param show_image whether to show the histology image in the background
#'
#' @returns A hex plot containing unique cell counts per spot
#' @export
#'
#' @examples
#' data("spatial_data_1")
#' deconv <- SpaceDeconv::deconvolute(spatial_data_1, method = "estimate")
#' SpaceDeconv::plot_cells_per_spot(deconv)
plot_cells_per_spot <- function(spatial_obj, plot_type = "spatial",
                                threshold = 0, spot_size = 1.5,
                                sample_id = "sample01", image_id = "lowres",
                                show_image = TRUE) {
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

  # remove unwanted columns
  res <- res[, !colnames(res) %in% c("in_tissue", "sample_id", "array_col", "array_row", "pxl_col_in_fullres", "pxl_row_in_fullres")]

  # count cells
  plot_data <- data.frame(spot = rownames(res), value = as.factor(apply(res, 1, sum)))

  df <- as.data.frame(cbind(SpatialExperiment::spatialCoords(spatial_obj), plot_data))

  return(make_baseplot(
    spe = spatial_obj, df = df, to_plot = "value",
    sample_id = sample_id, image_id = image_id,
    show_image = show_image, discrete = TRUE
  ))

  # if (plot_type == "bar") {
  #   plot <- ggplot2::ggplot(plot_data) +
  #     ggplot2::geom_bar(aes_string(x = "value"), stat = "count", fill = "black") +
  #     ggplot2::geom_vline(
  #       ggplot2::aes(xintercept = mean(unlist(plot_data["value"]))),
  #       color = "red",
  #       linetype = "dashed",
  #       size = 1
  #     ) +
  #     ggplot2::geom_hline(yintercept = 0, size = 1.2) +
  #     ggplot2::theme_classic() +
  #     ggplot2::theme(
  #       axis.title.x = ggplot2::element_blank(),
  #       axis.title.y = ggplot2::element_blank(),
  #       axis.line.y = ggplot2::element_blank(),
  #       axis.line.x = ggplot2::element_blank(),
  #     )
  # }
}

#' Function to plot deconvolution results
#'
#' Generate Hex Plot of a SpatialExperiment containing deconvolution results
#'
#' @param spe deconvolution result in Form of a SpatialExperiment
#' @param cell_type one or more celltype to plot
#' @param sample_id sample id to plot, default: "sample01"
#' @param image_id which image to plot, default: "lowres"
#' @param show_image logical, wether to display the image, default = TRUE
#' @param discrete logical, whether to scale the color discrete, default = FALSE
#'
#' @returns plot of cell type fractions
#'
#' @export
#' @examples
#' data("spatial_data_2")
#' deconv <- SpaceDeconv::deconvolute(spatial_data_2, method = "estimate")
#' SpaceDeconv::plot_celltype(deconv, cell_type = "estimate_immune.score")
plot_celltype <- function(spe, cell_type = NULL, sample_id = "sample01",
                          image_id = "lowres", show_image = TRUE, discrete = FALSE) {
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

  return(make_baseplot(spe, df,
    to_plot = cell_type, sample_id = sample_id,
    image_id = image_id, show_image = show_image,
    discrete = discrete
  ))

  # TODO
  # add facet wrap
  # add color selection
  # confirm the requested image is available in object
}


#' Function to plot deconvolution results
#'
#' Generate Hex Plot of a SpatialExperiment containing UMI counts
#'
#' @param spe deconvolution result in Form of a SpatialExperiment
#' @param sample_id sample id to plot, default: "sample01"
#' @param image_id which image to plot, default: "lowres"
#' @param show_image logical, wether to display the image, default = TRUE
#'
#' @returns plot of cell type fractions
#'
#' @export
#'
#' @examples
#' data("spatial_data_3")
#' deconv <- SpaceDeconv::deconvolute(spatial_data_3, method = "estimate")
#' plot_umi_count(deconv)
plot_umi_count <- function(spe, sample_id = "sample01", image_id = "lowres",
                           show_image = TRUE) {
  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }

  df <- as.data.frame(cbind(SpatialExperiment::spatialCoords(spe),
    nUMI = colSums(counts(spe))
  ))

  return(make_baseplot(spe,
    df,
    to_plot = "nUMI",
    sample_id = sample_id,
    image_id = image_id,
    show_image = show_image
  ))
}

###############
#### utils ####
###############


#' Render spatial hex plot
#'
#' @param spe SpatialExperiment with deconvolution results
#' @param df containing the annotation to be plotted
#' @param to_plot column of df to plot
#' @param sample_id sample of the SpatialExperiment to be plotted
#' @param image_id image id for background image
#' @param show_image whether to show the spatial image
#' @param discrete should the color scale be discrete? Defaut = FALSE
make_baseplot <- function(spe, df, to_plot, sample_id = "sample01",
                          image_id = "lowres", show_image = TRUE,
                          discrete = FALSE) {
  # scale coordinates with scalefactor
  df$pxl_col_in_fullres <- df$pxl_col_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sample_id, image_id = image_id)
  df$pxl_row_in_fullres <- df$pxl_row_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sample_id, image_id = image_id)

  # due to reasons, flip y axis by hand
  df$pxl_row_in_fullres <- df$pxl_row_in_fullres * -1

  # preparing the dataframe with sf, inserting points
  sf_points <- sf::st_as_sf(df, coords = c("pxl_col_in_fullres", "pxl_row_in_fullres"))

  # calculate spot distance
  spot_distance <- min(sqrt((df$pxl_col_in_fullres[1] - df$pxl_col_in_fullres[-1])^2 + (df$pxl_row_in_fullres[1] - df$pxl_row_in_fullres[-1])^2))

  # generate hexagons
  new_geom <- get_polygon_geometry(df, spot_distance)

  # no overwrite the points with hex polygons
  sf_poly <- sf::st_set_geometry(sf_points, new_geom)

  # extract image and dimensions
  img <- SpatialExperiment::imgRaster(spe, image_id = image_id)
  width <- dim(img)[2]
  height <- dim(img)[1]

  # initialize plot
  p <- ggplot()

  # add spatial image
  if (show_image) {
    p <- p + annotation_raster(img, xmin = 0, xmax = width, ymin = 0, ymax = -height)
  }

  # add hexagons
  p <- p +
    geom_sf(aes_string(fill = to_plot), lwd = 0, data = sf_poly) +
    coord_sf(xlim = c(0, width), ylim = c(0, -height)) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
    )

  # add color scale
  if (discrete) {
    p <- p + colorspace::scale_fill_discrete_sequential("Rocket")
  } else {
    p <- p + colorspace::scale_fill_continuous_sequential("Rocket")
  }

  return(p)
}

#' Build Hex Polygon Geometry
#'
#' @param x coordinate
#' @param y coordinate
#' @param dist distance of hexagons
get_hex_polygon <- function(x, y, dist) {
  angle <- seq(0, 2 * pi, length.out = 7)[-7] # angles of
  res <- cbind(
    x = x + sin(angle) * dist / 2,
    y = y + cos(angle) * dist / 2
  )
  return(rbind(res, res[1, ]))
}

# get coordinates for all hex polygons, move to utils
#' Build Full Hexagon Set for provided spots
#' @param grid Spot coordinates in df format, extracted from SpatialExperiment
#' @param dist distance between spots
get_polygon_geometry <- function(grid, dist) {
  res <- list()
  for (i in seq_len(nrow(grid))) {
    res <- c(
      res,
      list(sf::st_polygon(list(get_hex_polygon(
        grid$pxl_col_in_fullres[i],
        grid$pxl_row_in_fullres[i], dist
      ))))
    )
  }
  return(sf::st_sfc(res)) # Convert to 'simple feature collection'
}
