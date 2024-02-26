#' Ripley´s K Analysis for Cell Type Distribution (Beta)
#'
#' (BETA) Applies Ripley's K statistic to analyze the spatial distribution of a specified cell type within a `SpatialExperiment` object.
#' @param spe A `SpatialExperiment` object representing spatially resolved transcriptomics data.
#' @param cell_type The cell type of interest for which the Ripley's K statistic will be calculated.
#' @param method The deconvolution method used to determine the presence of the cell type within the spatial data.
#' @param threshold An optional cutoff value for determining cell type presence; if `NULL`, an antimode cutoff will be automatically determined based on the specified method.
#'
#' @return Returns the Ripley's K statistic for the specified cell type, providing insights into its spatial clustering within the dataset.
#'
#' @export
ripleys_k <- function(spe, cell_type, method, threshold = NULL) {
  coords <- spatialCoords(spe)

  # a <- antimode_cutoff(spe = spe, method = method, ), if threshold = NULL
  p <- presence(spe = spe, method = method, threshold = threshold)
  type <- as.factor(p[, cell_type])

  pp <- spatstat.geom::ppp(
    x = coords[, 2],
    y = coords[, 1],
    xrange = range(coords[, 2]),
    yrange = range(coords[, 1]),
    marks = type
  )
  k <- spatstat.explore::Kcross(pp, i = "TRUE", correction = "Ripley")
  return(k)
}


#' Largest radius and isotropic value
#'
#' @param k_functions a list with named ripley´s k statistics
#' @return largest radius and isotropic value in a list of k statistics
#'

get_largest_r_and_iso <- function(k_functions) {
  # Create reference values for r and iso
  largest_r <- -Inf
  largest_iso <- -Inf

  # Get the largest r and iso values
  for (i in seq_along(k_functions)) {
    k <- k_functions[[i]]
    max_r <- max(k$r)
    max_iso <- max(k$iso)

    if (max_r > largest_r) {
      largest_r <- max_r
    }
    if (max_iso > largest_iso) {
      largest_iso <- max_iso
    }
  }

  return(list(largest_r = largest_r, largest_iso = largest_iso))
}


## Use ggplot2 for the plotting function (Beta)

#' Combined Ripley´s K statistics in one plot (Beta)
#'
#' @param k_functions named list of k statistics from different SpatialExperiments
#' @returns combined Ripley´s K results in one plot
#' @export

plot_ripleys_k <- function(k_functions) {
  # Get largest r and iso value
  lims <- get_largest_r_and_iso(k_functions)

  # Convert the list to a data frame
  k_data <- data.frame(
    cell_type = names(k_functions),
    k_stats = I(k_functions),
    stringsAsFactors = FALSE
  )

  # Unnest the k_stats column
  k_data_unnested <- k_data %>%
    unnest(cols = k_stats)

  # Create an empty plot
  p <- ggplot() +
    ylim(0, lims$largest_iso) +
    xlim(0, lims$largest_r) +
    labs(x = "r", y = "K(r)", title = "Combined Ripley's K plot") +
    theme_minimal()

  # Poisson distribution
  p <- p + geom_line(data = k_data_unnested, aes(x = r, y = theo), color = "black", linetype = "dotted", size = 1)

  # Isotropic distribution of each cell type
  p <- p + geom_line(data = k_data_unnested, aes(x = r, y = iso, color = cell_type), size = 1)

  # Adjust the size of the text in the plot
  p <- p + theme(axis.text = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.title = element_text(size = 14), title = element_text(size = 14))

  return(p)
}
