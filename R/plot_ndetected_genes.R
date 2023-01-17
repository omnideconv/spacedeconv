#' Plot number of detected genes
#'
#' Generating a spatial plot of the number of detected genes of the SpatialExperiment
#'
#' @param spe deconvolution result in Form of a SpatialExperiment
#' @param palette colorspace palette (sequential)
#' @param transform_scale data transform_scaleation to use, "log"
#' @param reverse_palette reverse color palette
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
#' @param density whether to display a density distribution next to the spatial plot
#' @param save set TRUE to save plot
#' @param path specify directory to save plot, if NULL: saving at ~/spacedeconv
#' @param png_width when saving, png width in px
#' @param png_height when saving, png height in px
#'
#' @returns plot number of detected genes
#'
#' @export


plot_ndetected_genes <- function(spe, palette = "Mako", transform_scale = NULL,
                                 sample_id = "sample01", image_id = "lowres",
                                 reverse_palette = FALSE,
                                 show_image = FALSE, offset_rotation = FALSE,
                                 spot_size = 1, limits = NULL,
                                 smooth = FALSE, smoothing_factor = 1.5,
                                 title_size = 30, title = NULL, font_size = 20,
                                 legend_size = 40, density = TRUE,
                                 save = FALSE, path = NULL, png_width = 1500, png_height = 750) {
  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }

  gene_count <- counts(spe) >= 1

  df <- as.data.frame(cbind(SpatialExperiment::spatialCoords(spe),
    ndetected_genes = colSums(gene_count)
  ))

  return(make_baseplot(spe,
    df,
    to_plot = "ndetected_genes", sample_id = sample_id,
    image_id = image_id, show_image = show_image,
    offset_rotation = offset_rotation, palette = palette,
    transform_scale = transform_scale, reverse_palette = reverse_palette,
    spot_size = spot_size, limits = limits, title = title,
    smooth = smooth, smoothing_factor = smoothing_factor,
    title_size = title_size, font_size = font_size, legend_size = legend_size,
    density = density, save = save, path = path, png_width = png_width, png_height = png_height
  ))
}
