#' Quality control function
#'
#' Generates a html report containing quality metrics of the SpatialExperiment
#'
#' @param spe SpatialExperiment object
#' @param assay assay
#'
#' @returns html file with quality control metrics
#'
#' @export

qualitycontrol <- function(spe, assay = "counts") {
  # Add QC metrics
  spe <- addPerCellQC(spe)

  # number of total spots
  ncol(assay(spe, assay))

  # total number of genes
  dim(spe)[1]

  # Range of total UMIs per spot
  range(colData(spe)$sum)

  # % UMIs >= 500
  sum(colData(spe)$sum >= 500) / dim(spe)[2] * 100

  # % UMIs < 500
  sum(colData(spe)$sum < 500) / dim(spe)[2] * 100

  # Range of detected genes per spot
  range(colSums(assay(spe, assay) > 0))

  # Plot UMI
  plot_umi_count(spe, offset_rotation = TRUE)

  # Plot number of detected genes
  plot_ndetected_genes(spe, offset_rotation = TRUE)

  # Render html output
  rmarkdown::render(system.file("Qualitycontrol.rmd", package = "spacedeconv"), params = list(output_file = report.html))
}
