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

  gene_count <- counts(spe)
  for(i in 1:nrow(gene_count)){
    for(j in 1:ncol(gene_count)){
      if(gene_count[i,j] > 0){
        gene_count[i,j] <- 1
      }
    }
  }
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
