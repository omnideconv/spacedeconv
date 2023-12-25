#' Function to plot deconvolution results
#'
#' Generate Hex Plot of a SpatialExperiment containing deconvolution results
#'
#' @param spe deconvolution result in Form of a SpatialExperiment
#' @param cell_type one or more celltype to plot
#' @param palette colorspace palette (sequential)
#' @param transform_scale data transform_scaleation to use, "log"
#' @param reverse_palette reverse color palette
#' @param sample_id sample id to plot, default: "sample01"
#' @param image_id which image to plot, default: "lowres"
#' @param show_image logical, whether to display the image, default = TRUE
#' @param background background color
#' @param palette_type discrete, sequential or diverging
#' @param offset_rotation correct hex orientation for rotated visium image
#' @param spot_size increase (>1) or decrease (<1) the hex size
#' @param limits vector of color scale limits
#' @param smooth whether to smooth the plot
#' @param smoothing_factor kernel size factor (multiples of spot distance)
#' @param zoom zoom to the available spots
#' @param title_size font size of title
#' @param title set a custom title
#' @param font_size font size of legend
#' @param legend_size legend size in points
#' @param density whether to display a density distribution next to the spatial plot
#' @param save set TRUE to save plot
#' @param path specify directory to save plot, if NULL: saving at ~/spacedeconv
#' @param png_width when saving, png width in px
#' @param png_height when saving, png height in px
#' @param show_legend whether to show the legend
#'
#' @returns plot of cell type fractions
#'
#' @export
#' @examples
#' data("spatial_data_2")
#' deconv <- spacedeconv::deconvolute(spatial_data_2, method = "estimate")
#' spacedeconv::plot_celltype(deconv, cell_type = "estimate_immune.score")
plot_celltype <- function(spe, cell_type = NULL, palette = "Mako", transform_scale = NULL,
                          sample_id = "sample01", image_id = "lowres", reverse_palette = FALSE,
                          show_image = FALSE, background = NULL, palette_type = "sequential",
                          offset_rotation = FALSE, spot_size = 1, limits = NULL,
                          smooth = FALSE, smoothing_factor = 1.5, zoom = TRUE,
                          title_size = 30, title = NULL, font_size = 15, legend_size = 20, density = TRUE,
                          save = FALSE, path = NULL, png_width = 1500, png_height = 750,
                          show_legend = TRUE) {
  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if (is.null(cell_type)) {
    stop("Parameter 'cell_type' is null or missing, but is required")
  }

  # check that celltypes are present in object
  if (!all(cell_type %in% names(colData(spe))) && !cell_type %in% deconvolution_methods && !cell_type == "c2l" && !cell_type == "progeny" && !cell_type == "dorothea" && !cell_type == "cluster") {
    stop("Provides cell types are not present in SpatialExperiment")
  }

  spe <- filter_sample_id(spe, sample_id)

  df <- as.data.frame(cbind(SpatialExperiment::spatialCoords(spe), colData(spe)))


  # if a method is passed then make grid, otherwise, only one
  if (cell_type %in% deconvolution_methods || cell_type == "c2l" || cell_type == "progeny" || cell_type == "dorothea" || cell_type == "cluster") {
    plot <- make_baseplot(spe, df,
      palette = palette,
      to_plot = available_results(spe, method = cell_type)[1], sample_id = sample_id,
      image_id = image_id, show_image = show_image, background = background, zoom = zoom,
      palette_type = palette_type, offset_rotation = offset_rotation,
      transform_scale = transform_scale, reverse_palette = reverse_palette,
      spot_size = spot_size, limits = limits, title = title,
      smooth = smooth, smoothing_factor = smoothing_factor,
      title_size = title_size, font_size = font_size, legend_size = legend_size,
      density = density, save = save, path = path, png_width = png_width,
      png_height = png_height, show_legend = show_legend
    )

    for (result in available_results(spe, method = cell_type)[-1]) {
      plot <- plot + make_baseplot(spe, df,
        palette = palette,
        to_plot = result, sample_id = sample_id,
        image_id = image_id, show_image = show_image, background = background, zoom = zoom,
        palette_type = palette_type, offset_rotation = offset_rotation,
        transform_scale = transform_scale, reverse_palette = reverse_palette,
        spot_size = spot_size, limits = limits, title = title,
        smooth = smooth, smoothing_factor = smoothing_factor,
        title_size = title_size, font_size = font_size, legend_size = legend_size,
        density = density, save = save, path = path, png_width = png_width,
        png_height = png_height, show_legend = show_legend
      )
    }

    return(plot)
  } else {
    return(make_baseplot(spe, df,
      palette = palette,
      to_plot = cell_type, sample_id = sample_id,
      image_id = image_id, show_image = show_image, background = background, zoom = zoom,
      palette_type = palette_type, offset_rotation = offset_rotation,
      transform_scale = transform_scale, reverse_palette = reverse_palette,
      spot_size = spot_size, limits = limits, title = title,
      smooth = smooth, smoothing_factor = smoothing_factor,
      title_size = title_size, font_size = font_size, legend_size = legend_size,
      density = density, save = save, path = path, png_width = png_width,
      png_height = png_height, show_legend = show_legend
    ))
  }


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
#' @param palette colorspace palette (sequential)
#' @param transform_scale data transform_scaleation to use, "log"
#' @param reverse_palette reverse color palette
#' @param sample_id sample id to plot, default: "sample01"
#' @param image_id which image to plot, default: "lowres"
#' @param show_image logical, wether to display the image, default = TRUE
#' @param background custom background color
#' @param zoom zoom to the available spots
#' @param offset_rotation correct hex orientation for rotated visium image
#' @param spot_size increase (>1) or decrease (<1) the hex size
#' @param limits vector of color scale limits
#' @param smooth whether to smooth the plot
#' @param smoothing_factor kernel size factor (multiples of spot distance)
#' @param title_size font size of title
#' @param title set a custom title
#' @param font_size font size of legend
#' @param legend_size legend size in points
#' @param density whether to display a density distribution next to the spatial plot
#' @param save set TRUE to save plot
#' @param path specify directory to save plot, if NULL: saving at ~/spacedeconv
#' @param png_width when saving, png width in px
#' @param png_height when saving, png height in px
#' @param show_legend whether to show the legend
#'
#'
#' @returns plot of cell type fractions
#'
#' @export
#'
#' @examples
#' data("spatial_data_3")
#' deconv <- spacedeconv::deconvolute(spatial_data_3, method = "estimate")
#' plot_umi_count(deconv)
plot_umi_count <- function(spe, palette = "Mako", transform_scale = NULL,
                           sample_id = "sample01", image_id = "lowres",
                           reverse_palette = FALSE, zoom = TRUE,
                           show_image = FALSE, background = NULL, offset_rotation = FALSE,
                           spot_size = 1, limits = NULL,
                           smooth = FALSE, smoothing_factor = 1.5,
                           title_size = 30, title = NULL, font_size = 15,
                           legend_size = 20, density = TRUE,
                           save = FALSE, path = NULL, png_width = 1500, png_height = 750,
                           show_legend = TRUE) {
  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }

  spe <- filter_sample_id(spe, sample_id)

  df <- as.data.frame(cbind(SpatialExperiment::spatialCoords(spe),
    nUMI = colSums(counts(spe))
  ))

  return(make_baseplot(spe,
    df,
    to_plot = "nUMI", sample_id = sample_id,
    image_id = image_id, show_image = show_image, background = background, zoom = zoom,
    offset_rotation = offset_rotation, palette = palette,
    transform_scale = transform_scale, reverse_palette = reverse_palette,
    spot_size = spot_size, limits = limits, title = title,
    smooth = smooth, smoothing_factor = smoothing_factor,
    title_size = title_size, font_size = font_size, legend_size = legend_size,
    density = density, save = save, path = path, png_width = png_width, png_height = png_height,
    show_legend = show_legend
  ))
}

#' Function to plot deconvolution results
#'
#' Generate Hex Plot of a SpatialExperiment containing the most abundant cell types
#'
#' @param spe deconvolution result in Form of a SpatialExperiment
#' @param method select which results should be displayed
#' @param cell_type one or more celltype to plot, NULL for all
#' @param remove vector of cell types to be removed from the plot
#' @param min_spot minimum number of spots the cell-type has to be present in
#' @param palette colorspace palette (sequential)
# #' @param transform_scale data transform_scaleation to use, "log"
#' @param reverse_palette reverse color palette
#' @param sample_id sample id to plot, default: "sample01"
#' @param image_id which image to plot, default: "lowres"
#' @param show_image logical, whether to display the image, default = TRUE
#' @param background custom background color
#' @param zoom zoom to the available spots
#' @param palette_type logical, whether to scale the color palette_type, default = FALSE
#' @param offset_rotation correct hex orientation for rotated visium image
#' @param spot_size increase (>1) or decrease (<1) the hex size
# #' @param limits vector of color scale limits
# #' @param smooth whether to smooth the plot
# #' @param smoothing_factor kernel size factor (multiples of spot distance)
#' @param title_size font size of title
#' @param title set a custom title
#' @param font_size font size of legend
#' @param legend_size legend size in points
#' @param density whether to display a density distribution next to the spatial plot
#' @param save set TRUE to save plot
#' @param path specify directory to save plot, if NULL: saving at ~/spacedeconv
#' @param png_width when saving, png width in px
#' @param png_height when saving, png height in px
#' @param show_legend whether to show the legend
#' @param min_abundance minimum abundance of celltypes to be included in the analysis
#'
#' @returns plot of cell type fractions
#'
#' @export
plot_most_abundant <- function(spe, method = NULL, cell_type = NULL, remove = NULL, min_spot = 0, palette = "Mako", # transform_scale = NULL,
                               sample_id = "sample01", image_id = "lowres", reverse_palette = FALSE,
                               show_image = FALSE, background = NULL, zoom = TRUE, palette_type = "discrete",
                               offset_rotation = FALSE, spot_size = 1, # limits = NULL,
                               # smooth = FALSE, smoothing_factor = 1.5,
                               title_size = 30, font_size = 15, legend_size = 20,
                               density = FALSE, save = FALSE, path = NULL,
                               png_width = 1500, png_height = 750, title = NULL,
                               show_legend = TRUE, min_abundance = 0.01) {
  # checks
  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if (is.null(method)) {
    stop("Parameter 'method' is null or missing, but is required")
  }

  if (!is.null(method)) {
    available <- available_results(spe)[startsWith(available_results(spe), method)]
  } else {
    available <- available_results(spe)
  }

  if (!is.null(cell_type)) {
    available <- cell_type
  }

  if (!is.null(remove)) {
    available <- available[!available %in% remove]
  }

  # filter sample from object
  spe <- filter_sample_id(spe, sample_id)

  # create df with data to plot
  df <- as.data.frame(colData(spe))[, available, drop = FALSE]
  df <- df[, !names(df) %in% c("in_tissue", "array_row", "array_col", "sample_id"), drop = FALSE]

  # remove all columns not numeric
  df <- df[, unlist(lapply(df, is.numeric)), drop = FALSE]

  # handle min_abundance: set all other to 0
  df[df < min_abundance] <- 0

  # ensure min_spot parameter
  if (min_spot > 0) {
    # compute how many spots have >0 values for each celltype
    n_above_zero <- sapply(df, function(column) {
      sum(column != 0)
    })

    # subset df again if celltypes to sparse
    df <- df[, names(n_above_zero)[n_above_zero > min_spot]]
  }


  # list of the mostAbundant celltype per spot
  res <- colnames(df)[max.col(df)]

  # update the rows with all zero rows to specific string
  res[rowSums(df) == 0] <- "Not enough Data"

  # append result to df
  df2 <- as.data.frame(cbind(SpatialExperiment::spatialCoords(spe), mostAbundant = res))

  df2$pxl_col_in_fullres <- as.numeric(df2$pxl_col_in_fullres)
  df2$pxl_row_in_fullres <- as.numeric(df2$pxl_row_in_fullres)


  return(make_baseplot(
    spe = spe, df = df2, to_plot = "mostAbundant", palette = palette,
    sample_id = sample_id, image_id = image_id, background = background, zoom = zoom,
    reverse_palette = reverse_palette, show_image = show_image,
    offset_rotation = offset_rotation, spot_size = spot_size,
    title_size = title_size, palette_type = palette_type,
    font_size = font_size, legend_size = legend_size, density = density,
    save = save, path = path, png_width = png_width, png_height = png_height,
    title = title, show_legend = show_legend,
  ))
}

#' Plot celltype presence absence
#'
#' @param spe deconvolution result in Form of a SpatialExperiment
#' @param cell_type celltype to plot
#' @param threshold fraction threshold, if NULL: calculated internally
#' @param palette colorspace palette (sequential)
#' @param transform_scale data transform_scaleation to use, "log"
#' @param reverse_palette reverse color palette
#' @param background custom background color
#' @param zoom zoom to the available spots
#' @param sample_id sample id to plot, default: "sample01"
#' @param image_id which image to plot, default: "lowres"
#' @param show_image logical, wether to display the image, default = TRUE
#' @param offset_rotation correct hex orientation for rotated visium image
#' @param spot_size increase (>1) or decrease (<1) the hex size
#' @param limits vector of color scale limits
#' @param smooth whether to smooth the plot
#' @param smoothing_factor kernel size factor (multiples of spot distance)
#' @param title_size font size of title
#' @param title set a custom title
#' @param font_size font size of legend
#' @param legend_size legend size in points
#' @param save set TRUE to save plot
#' @param path specify directory to save plot, if NULL: saving at ~/spacedeconv
#' @param png_width when saving, png width in px
#' @param png_height when saving, png height in px
#' @param show_legend whether to show the legend
#'
#' @returns plot of a celltypes presence/absence using a threshold
#'
#' @export
plot_celltype_presence <- function(spe, cell_type = NULL, threshold = NULL,
                                   palette = "Mako", transform_scale = NULL,
                                   sample_id = "sample01", image_id = "lowres",
                                   reverse_palette = FALSE, background = NULL, zoom = TRUE,
                                   show_image = FALSE, offset_rotation = FALSE,
                                   spot_size = 1, limits = NULL,
                                   smooth = FALSE, smoothing_factor = 1.5,
                                   title_size = 30, title = NULL, font_size = 15,
                                   legend_size = 20, save = FALSE, path = NULL,
                                   png_width = 1500, png_height = 750, show_legend = TRUE) {
  spe <- filter_sample_id(spe, sample_id)

  df <- as.data.frame(cbind(SpatialExperiment::spatialCoords(spe), colData(spe)))

  # extract method from celltype
  method <- unlist(strsplit(cell_type, "_"))[1]

  if (is.null(threshold)) {
    threshold <- antimode_cutoff(spe, method)[cell_type]
    message("Calculated threshold for ", cell_type, ": ", threshold)
  }

  # calculate presence
  presence <- presence(spe, method, threshold)[, cell_type]
  df <- cbind(df, presence = presence)

  if (is.null(title)) {
    title <- paste0("presence_", cell_type)
  }

  return(make_baseplot(
    spe = spe, df = df, to_plot = "presence", palette = palette,
    transform_scale = transform_scale, sample_id = sample_id,
    image_id = image_id, reverse_palette = reverse_palette,
    show_image = show_image, offset_rotation = offset_rotation,
    spot_size = spot_size, limits = limits, smooth = smooth,
    smoothing_factor = smoothing_factor, title_size = title_size,
    font_size = font_size, legend_size = legend_size, background = background, zoom = zoom,
    density = FALSE, palette_type = "discrete", save = save, path = path,
    png_width = png_width, png_height = png_height, title = title, show_legend = show_legend
  ))
}

#' Plot celltype fraction comparison
#'
#' @param spe deconvolution result in Form of a SpatialExperiment
#' @param cell_type_1 celltype to plot
#' @param cell_type_2 celltype to plot
#' @param palette colorspace palette (sequential)
#' @param transform_scale data transform_scaleation to use, "log"
#' @param reverse_palette reverse color palette
#' @param sample_id sample id to plot, default: "sample01"
#' @param image_id which image to plot, default: "lowres"
#' @param background custom background color
#' @param zoom zoom to the available spots
#' @param show_image logical, wether to display the image, default = TRUE
#' @param offset_rotation correct hex orientation for rotated visium image
#' @param spot_size increase (>1) or decrease (<1) the hex size
#' @param limits vector of color scale limits
#' @param smooth whether to smooth the plot
#' @param smoothing_factor kernel size factor (multiples of spot distance)
#' @param title_size font size of title
#' @param title set a custom title
#' @param font_size font size of legend
#' @param legend_size legend size in points
#' @param density whether to display a density distribution next to the spatial plot
#' @param palette_type "discrete", "sequenatial", "diverging"
#' @param save set TRUE to save plot
#' @param path specify directory to save plot, if NULL: saving at ~/spacedeconv
#' @param png_width when saving, png width in px
#' @param png_height when saving, png height in px
#' @param show_legend whether to show the legend
#'
#' @returns plot of a celltypes presence/absence using a threshold
#'
#' @export
plot_comparison <- function(spe, cell_type_1 = NULL, cell_type_2 = NULL,
                            palette = "Blue-Red", transform_scale = NULL,
                            sample_id = "sample01", image_id = "lowres",
                            reverse_palette = FALSE, background = NULL, zoom = TRUE,
                            show_image = FALSE, offset_rotation = FALSE,
                            spot_size = 1, limits = NULL,
                            smooth = FALSE, smoothing_factor = 1.5,
                            title_size = 30, title = NULL, font_size = 15,
                            legend_size = 20, palette_type = "diverging", density = TRUE,
                            save = FALSE, path = NULL, png_width = 1500, png_height = 750,
                            show_legend = TRUE) {
  spe <- filter_sample_id(spe, sample_id)

  df <- as.data.frame(cbind(SpatialExperiment::spatialCoords(spe), colData(spe)))

  # comparison <-df[, cell_type_1] - df[, cell_type_2]
  #
  # cmean <- mean(comparison)
  # csd <- sd(comparison)
  #
  # zcomparison <- (comparison-cmean)/csd
  # comparison <- zcomparison


  comparison <- (df[, cell_type_1] + 1) / (df[, cell_type_2] + 1)
  # comparison <- comparison - 1
  comparison <- log(comparison)
  comparison[is.infinite(comparison)] <- NA # ?

  df <- cbind(df, comparison = comparison)

  # custom title
  if (is.null(title)) {
    title <- paste0("comparison", cell_type_1, "_", "cell_type_2")
  }

  return(make_baseplot(
    spe = spe, df = df, to_plot = "comparison", palette = palette,
    transform_scale = transform_scale, sample_id = sample_id,
    image_id = image_id, reverse_palette = reverse_palette,
    show_image = show_image, offset_rotation = offset_rotation,
    spot_size = spot_size, limits = limits, smooth = smooth,
    smoothing_factor = smoothing_factor, title_size = title_size,
    font_size = font_size, legend_size = legend_size, background = background, zoom = zoom,
    density = density, palette_type = palette_type, save = save, path = path,
    png_width = png_width, png_height = png_height, title = title, show_legend = show_legend
  ))
}


#' Function to plot gene expression
#'
#' Generate Hex Plot of a SpatialExperiment containing deconvolution results
#'
#' @param spe deconvolution result in Form of a SpatialExperiment
#' @param gene gene to plot
#' @param assay assay to extract gene expression, default = "counts"
#' @param palette colorspace palette (sequential)
#' @param transform_scale data transform_scaleation to use, "log"
#' @param reverse_palette reverse color palette
#' @param sample_id sample id to plot, default: "sample01"
#' @param image_id which image to plot, default: "lowres"
#' @param show_image logical, whether to display the image, default = TRUE
#' @param background custom background color
#' @param zoom zoom to the available spots
#' @param palette_type discrete, sequential or diverging
#' @param offset_rotation correct hex orientation for rotated visium image
#' @param spot_size increase (>1) or decrease (<1) the hex size
#' @param limits vector of color scale limits
#' @param smooth whether to smooth the plot
#' @param smoothing_factor kernel size factor (multiples of spot distance)
#' @param show_legend whether to show the legend
#' @param title_size font size of title
#' @param title set a custom title
#' @param font_size font size of legend
#' @param legend_size legend size in points
#' @param density whether to display a density distribution next to the spatial plot
#' @param save set TRUE to save plot
#' @param path specify directory to save plot, if NULL: saving at ~/spacedeconv
#' @param png_width when saving, png width in px
#' @param png_height when saving, png height in px
#'
#' @returns plot of cell type fractions
#'
#' @export
#' @examples
#' data("spatial_data_2")
#' deconv <- spacedeconv::deconvolute(spatial_data_2, method = "estimate")
#' spacedeconv::plot_celltype(deconv, cell_type = "estimate_immune.score")
plot_gene <- function(spe, gene = NULL, assay = "counts", palette = "Mako", transform_scale = NULL,
                      sample_id = "sample01", image_id = "lowres", reverse_palette = FALSE,
                      show_image = FALSE, background = NULL, zoom = TRUE, palette_type = NULL, # sequential
                      offset_rotation = FALSE, spot_size = 1, limits = NULL,
                      smooth = FALSE, smoothing_factor = 1.5, show_legend = TRUE,
                      title_size = 30, title = NULL, font_size = 15, legend_size = 20, density = TRUE,
                      save = FALSE, path = NULL, png_width = 1500, png_height = 750) {
  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if (is.null(gene)) {
    stop("Parameter 'gene' is null or missing, but is required")
  }

  # check that celltypes are present in object
  if (!gene %in% rownames(spe)) {
    stop("Provides gene is not present in SpatialExperiment")
  }

  spe <- filter_sample_id(spe, sample_id)

  df <- as.data.frame(cbind(SpatialExperiment::spatialCoords(spe), gene = SummarizedExperiment::assay(spe, assay)[gene, ]))

  # set gene name as title
  if (is.null(title)) {
    title <- gene
  }

  return(make_baseplot(spe, df,
    title = title,
    palette = palette,
    to_plot = "gene", sample_id = sample_id,
    image_id = image_id, show_image = show_image,
    palette_type = palette_type, offset_rotation = offset_rotation,
    transform_scale = transform_scale, reverse_palette = reverse_palette,
    spot_size = spot_size, limits = limits, background = background, zoom = zoom,
    smooth = smooth, smoothing_factor = smoothing_factor, show_legend = show_legend,
    title_size = title_size, font_size = font_size, legend_size = legend_size,
    density = density, save = save, path = path, png_width = png_width, png_height = png_height
  ))

  # TODO
  # add facet wrap
  # add color selection
  # confirm the requested image is available in object
}


###############
#### utils ####
###############

#' Render spatial hex plot
#'
#' @param spe SpatialExperiment with deconvolution results
#' @param df containing the annotation to be plotted
#' @param to_plot column of df to plot
#' @param palette colorspace palette (sequential)
#' @param transform_scale data transform_scaleation to use, "log"
#' @param reverse_palette reverse color palette
#' @param sample_id sample of the SpatialExperiment to be plotted
#' @param image_id image id for background image
#' @param show_image whether to show the spatial image
#' @param background custom background color
#' @param zoom zoom to the available spots
#' @param palette_type "discrete", "sequential", "diverging"
#' @param offset_rotation correct hex orientation for rotated visium image
#' @param spot_size increase (>1) or decrease (<1) the hex size
#' @param limits vector of color scale limits
#' @param smooth whether to smooth the plot
#' @param smoothing_factor kernel size factor (multiples of spot distance)
#' @param title_size font size of title
#' @param title set a custom title
#' @param font_size font size of legend
#' @param legend_size legend size in points
#' @param density whether to display a density distribution next to the spatial plot
#' @param save set TRUE to save plot
#' @param path specify directory to save plot, if NULL: saving at ~/spacedeconv
#' @param png_width when saving, png width in px
#' @param png_height when saving, png height in px
#' @param show_legend whether to show the legend
make_baseplot <- function(spe, df, to_plot, palette = "Mako", transform_scale = NULL,
                          sample_id = "sample01", reverse_palette = FALSE,
                          image_id = "lowres", show_image = FALSE, background = NULL, zoom = TRUE,
                          palette_type = NULL, offset_rotation = FALSE, spot_size = 1,
                          limits = NULL, smooth = FALSE, smoothing_factor = 1.5,
                          title_size = 30, title = NULL, font_size = 15, legend_size = 20, density = TRUE,
                          save = FALSE, path = NULL, png_width = 1500, png_height = 750, show_legend = TRUE) {
  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if (is.null(df)) {
    stop("Parameter 'df' is null or missing, but is required")
  }

  if (is.null(to_plot)) {
    stop("Parameter 'to_plot' is null or missing, but is required")
  }

  # scale coordinates with scalefactor
  df$pxl_col_in_fullres <- df$pxl_col_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sample_id, image_id = image_id)
  df$pxl_row_in_fullres <- df$pxl_row_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sample_id, image_id = image_id)

  # due to reasons, flip y axis by hand
  df$pxl_row_in_fullres <- df$pxl_row_in_fullres * -1

  # calculate scaling factor
  scaling_offset <- 1.165 # 1.154701 # 1/cos((30/360)*2*pi)

  # calculate spot distance
  spot_distance <- min(sqrt((df$pxl_col_in_fullres[1] - df$pxl_col_in_fullres[-1])^2 + (df$pxl_row_in_fullres[1] - df$pxl_row_in_fullres[-1])^2)) * spot_size * scaling_offset

  # smooth if requested
  if (smooth) {
    df[[to_plot]] <- smooth_celltype(df, spot_distance = spot_distance, smoothing_factor = smoothing_factor, cell_type = to_plot)
  }

  # tranform scale, if NULL is provided set to none for the switch to work
  transform_scale <- tolower(ifelse(is.null(transform_scale), "none", transform_scale))
  df[[to_plot]] <- switch(transform_scale,
    "ln" = log((df[[to_plot]] - min(df[[to_plot]])) + 1),
    "log10" = log10((df[[to_plot]] - min(df[[to_plot]])) + 1),
    "log2" = log2((df[[to_plot]] - min(df[[to_plot]])) + 1),
    "sqrt" = sqrt(df[[to_plot]]),
    "log" = log((df[[to_plot]] - min(df[[to_plot]])) + 1),
    df[[to_plot]]
  ) # Default case: no transformation
  transform_suffix <- if (transform_scale %in% c("ln", "log10", "log2", "sqrt", "log")) transform_scale else ""


  # manually fix limits, overwrite values to prevent NA
  if (!is.null(limits)) {
    df[df[, to_plot] < limits[1], to_plot] <- limits[1] # lower limit
    df[df[, to_plot] > limits[2], to_plot] <- limits[2] # upper limit
  }



  # Check if plot is smoothed
  if (smooth) {
    smooth_suffix <- "smoothed"
  } else {
    smooth_suffix <- ""
  }

  # handle legend title
  if (is.null(title)) {
    # Construct legend_title based on transformation and smoothing
    legend_title <- as.character(to_plot)
    if (transform_suffix != "") {
      legend_title <- paste0(legend_title, "_", transform_suffix)
    }
    if (smooth_suffix != "") {
      legend_title <- paste0(legend_title, "_", smooth_suffix)
    }
  } else {
    # custom title overwrites everything
    legend_title <- title
  }

  # preparing the dataframe with sf, inserting points
  sf_points <- sf::st_as_sf(df, coords = c("pxl_col_in_fullres", "pxl_row_in_fullres"))

  # generate hexagons
  new_geom <- get_polygon_geometry(df, spot_distance, offset_rotation = offset_rotation)

  # now overwrite the points with hex polygons
  sf_poly <- sf::st_set_geometry(sf_points, new_geom)

  # Determine the palette type based on the palette name
  if (is.null(palette_type)) {
    palette_type <- get_palette_type(palette)
  }

  # check discrete and, if yes, remove the hexagons which should not be plotted

  if (palette_type == "discrete") {
    tmp <- as.data.frame(sf_poly)
    if (is.logical(tmp[, to_plot])) {
      sf_poly <- sf_poly[tmp[, to_plot], ]
    }
  }

  # extract image and dimensions
  img <- SpatialExperiment::imgRaster(spe, image_id = image_id)
  width <- dim(img)[2]
  height <- dim(img)[1]

  # initialize plot
  p <- ggplot()

  # add spatial image
  if (show_image) {
    p <- p + annotation_raster(img, xmin = 0, xmax = width, ymin = 0, ymax = -height)
    # p <- p + annotation_raster(img, xmin = min(df$pxl_col_in_fullres), xmax = max(df$pxl_col_in_fullres), ymin = max(df$pxl_row_in_fullres), ymax = min(df$pxl_row_in_fullres))
  }

  # add hexagons
  p <- p +
    geom_sf(aes_string(fill = to_plot), lwd = 0, color = NA, data = sf_poly) +
    # coord_sf(xlim = c(0, width), ylim = c(0, -height)) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = title_size, hjust = 0.5),
      legend.text = element_text(size = font_size),
      legend.key.size = unit(legend_size, "points")
    ) +
    ggplot2::labs(title = legend_title, fill = element_blank())

  # zoom if requested
  if (zoom) {
    p <- p + coord_sf(xlim = c(min(df$pxl_col_in_fullres), max(df$pxl_col_in_fullres)), ylim = c(max(df$pxl_row_in_fullres), min(df$pxl_row_in_fullres)))
    # p <- p + coord_sf(xlim = c(min_col, max_col), ylim = c(max_row, min_row))
  } else {
    p <- p + coord_sf(xlim = c(0, width), ylim = c(0, -height))
  }

  # show legend?
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  # add custom background
  if (!is.null(background)) {
    # check that color actually is a color
    p <- p + theme(panel.background = element_rect(fill = background))
  } else {
    p <- p + theme(panel.background = element_blank())
  }

  # choose the palette type/coloring
  # First manual, then R Color Brewer, then colorspace
  if (is.vector(palette) && all(sapply(palette, function(x) is.character(x) && colorspace::is_hexcolor(x)))) {
    # manual palette
    p <- p + ggplot2::scale_fill_manual(values = palette)
  } else if (is.character(palette)) {
    # RColorBrewer
    if (palette %in% rownames(RColorBrewer::brewer.pal.info)) {
      # Number of unique values to be plotted
      num_values <- length(unique(df[[to_plot]]))

      max_colors <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]

      # if not enough colors then interpolate
      if (num_values > max_colors) {
        cli::cli_alert_info("Palette too small, interpolating colors!")

        # Interpolate the palette to get enough colors
        brewer_palette <- RColorBrewer::brewer.pal(max_colors, palette)
        interpolated_palette <- colorRampPalette(brewer_palette)
        palette_function <- interpolated_palette(num_values)
        p <- p + ggplot2::scale_fill_manual(values = palette_function)
      } else {
        p <- p + scale_fill_brewer(palette = palette) # reverse!
      }
    } else {
      # Colorspace

      if (is.factor(df[[to_plot]]) || is.character(df[[to_plot]]) || is.logical(df[[to_plot]]) || palette_type == "discrete") {
        # p <- p + colorspace::scale_fill_discrete_sequential("Inferno", rev = reverse_palette, limits = limits)
        # manual fix !!!
        if (palette_type == "sequential" || palette_type == "discrete") {
          pal <- function(n) {
            colorspace::sequential_hcl(n, palette, rev = reverse_palette)
          }
        } else if (palette_type == "diverging") {
          pal <- function(n) {
            colorspace::diverging_hcl(n, palette, rev = reverse_palette)
          }
        } else if (palette_type == "qualitative") {
          pal <- function(n) {
            colorspace::qualitative_hcl(n, palette, rev = reverse_palette)
          }
        } else {
          print("Error while adding color palette")
        }

        p <- p + ggplot2::discrete_scale(aesthetics = "fill", "manual", pal)
      } else if (palette_type == "sequential") {
        p <- p + colorspace::scale_fill_continuous_sequential(palette, rev = reverse_palette, limits = limits)
      } else if (palette_type == "diverging") {
        p <- p + colorspace::scale_fill_continuous_diverging(palette, rev = reverse_palette, limits = limits)
      }
    }
  } else {
    stop("Invalid palette input. It must be either a named palette or a vector of color hex codes.")
  }

  # create density plot if requested
  suppressMessages(
    if (density && palette_type != "discrete") {
      data <- data.frame(values = sf_poly[[to_plot]], id = rep(to_plot, nrow(sf_poly)))
      density <- ggplot2::ggplot(data, mapping = ggplot2::aes_string(x = "values", y = "id")) + # fill... see ggridges docs
        # ggplot2::geom_density() +
        ggridges::geom_density_ridges() + # _gradient()
        # ggplot2::scale_fill_viridis_c() +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
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
          axis.title = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(size = 14)
        )
      # ggplot2::ylim(0, max(data["values"]))

      # cowplot::plot_grid(spatial, density)
      plot <- ggpubr::ggarrange(p, density, ncol = 2, widths = c(2, 1)) # add functions to pkg.R
      # plot <- grid::grid.draw(plot) # add functions to pkg.R
    } else {
      plot <- p
    }
  )

  if (save) {
    # check if provided path works
    if (is.null(path) || !file.exists(path)) {
      # set default
      if (!file.exists("~/spacedeconvResults")) {
        dir.create("~/spacedeconvResults")
      }
      path <- normalizePath("~/spacedeconvResults")
    }
    save_plot(plot, to_plot, path, png_width, png_height)
  }


  plot # did not work with return ()
}

#' Build Hex Polygon Geometry
#'
#' @param x coordinate
#' @param y coordinate
#' @param dist distance of hexagons
#' @param offset_rotation correct hex orientation for rotated visium image
get_hex_polygon <- function(x, y, dist, offset_rotation = FALSE) {
  # offset rotation for visium image
  offset <- 0
  if (offset_rotation) {
    offset <- 0.5 # rotate by 30 degrees
  }

  angle <- seq(0, 2 * pi, length.out = 7)[-7] # angles of
  res <- cbind(
    x = x + sin(angle + offset) * dist / 2,
    y = y + cos(angle + offset) * dist / 2
  )
  return(rbind(res, res[1, ]))
}

# get coordinates for all hex polygons, move to utils
#' Build Full Hexagon Set for provided spots
#' @param grid Spot coordinates in df format, extracted from SpatialExperiment
#' @param dist distance between spots
#' @param offset_rotation correct hex orientation for rotated visium image
get_polygon_geometry <- function(grid, dist, offset_rotation = FALSE) {
  res <- list()
  for (i in seq_len(nrow(grid))) {
    res <- c(
      res,
      list(sf::st_polygon(list(get_hex_polygon(
        grid$pxl_col_in_fullres[i],
        grid$pxl_row_in_fullres[i], dist,
        offset_rotation = offset_rotation
      ))))
    )
  }
  return(sf::st_sfc(res)) # Convert to 'simple feature collection'
}
#' Smooth spot annotation
#'
#' @param df DataFrame containing spot coordinates and the value to be smoothed
#' @param spot_distance distance between two spots
#' @param smoothing_factor multiplied with spot distance to find close spots
#' @param cell_type value to be smoothed
smooth_celltype <- function(df, spot_distance, smoothing_factor = 1.5, cell_type = NULL) {
  new_values <- vector()

  # for all spots get the spots in distance and calculate mean value
  for (spot in rownames(df)) {
    point <- df[spot, ][c("pxl_col_in_fullres", "pxl_row_in_fullres")]

    spots_in_distance <- sqrt((point$pxl_col_in_fullres - df$pxl_col_in_fullres)^2 + (point$pxl_row_in_fullres - df$pxl_row_in_fullres)^2)
    names(spots_in_distance) <- rownames(df)
    spots_in_distance <- spots_in_distance[spots_in_distance <= spot_distance * smoothing_factor]

    new_values <- c(new_values, mean(df[names(spots_in_distance), cell_type], na.rm = TRUE))
  }

  return(new_values)
}


#' Save Plot to path
#'
#' @param plot ggplot
#' @param to_plot celltype name
#' @param path where the plot should be stored
#' @param png_width when saving, png width in px
#' @param png_height when saving, png height in px
#'
save_plot <- function(plot, to_plot, path, png_width, png_height) {
  filename <- paste0(path, "/", to_plot, ".png")
  png(width = png_width, height = png_height, units = "px", filename = filename)
  grid::grid.draw(plot)
  dev.off()
}

#' Filter SPE to contain only one sample ID
#' @param spe SpatialExperiment
#' @param sample_id sample_id
#' @export
filter_sample_id <- function(spe, sample_id) {
  if (is.null(spe)) {
    cli::cli_alert_danger("Spatial Object not provided")
    stop()
  }

  # check if sample id is provided, if no and only one available use this one
  if (is.null(sample_id)) {
    if (length(unique(spe$sample_id)) == 1) {
      cli::cli_alert_info("No sample ID provided, using the only one available")
      sample_id <- unique(spe$sample_id)
    } else {
      cli::cli_alert_danger("Multiple Sample Ids in Objet, please select one")
      stop()
    }
  }

  # remove columns not in this sample
  spe <- spe[, spe$sample_id == sample_id]

  return(spe)
}

#' Determine the Type of a Given Color Palette
#'
#' This function identifies the type of a color palette, supporting palettes
#' from both the `colorspace` and `RColorBrewer` packages. It is used to
#' streamline the process of palette type identification for use in plotting functions.
#'
#' @param palette A character string specifying the name of the color palette.
#'                The name should be a valid palette name from either the
#'                `colorspace` or `RColorBrewer` package.
#'
#' @return A character string representing the type of the palette.
#'         Possible return values are "sequential", "diverging", or "qualitative".
#'         The function will stop and display an error message if an unknown palette
#'         name is provided.
#'
get_palette_type <- function(palette) {
  palettes_colorspace <- colorspace::hcl_palettes()
  if (palette %in% rownames(palettes_colorspace)) {
    palette_info <- palettes_colorspace[palette, ]
    return(tolower(palette_info$type))
  } else if (palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    palette_info <- RColorBrewer::brewer.pal.info[palette, ]
    brewer_type <- tolower(palette_info$category)
    type_mapping <- c("div" = "diverging", "seq" = "sequential", "qual" = "qualitative")
    return(type_mapping[brewer_type])
  } else {
    stop("Unknown palette. Please use a palette from colorspace or RColorBrewer.")
  }
}
