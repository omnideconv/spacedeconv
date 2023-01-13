## Quality control function

qualitycontrol <- function(spe) {
  # Add QC metrics
  spe <- addPerCellQC(spe)

  # number of total spots
  ncol(counts(spe))

  # total number of genes
  dim(spe)[1]

  # Range of total UMIs per spot
  range(colData(spe)$sum)

  # % UMIs >= 500
  sum(colData(spe)$sum >= 500) / dim(spe)[1]

  # % UMIs < 500
  sum(colData(spe)$sum < 500) / dim(spe)[1]

  # Range of detected genes per spot
  range(colSums(counts(spe) > 0))

  # Plot UMI
  plot_umi_count(spe, offset_rotation = T)

  # Plot number of detected genes
  plot_ndetected_genes(spe, offset_rotation = T)

  # Render html output
  rmarkdown::render('./Qualitycontrol.rmd', params = list(output_file = report.html))
  #rmarkdown::render(system.file("Qualitycontrol.rmd", package = "spacedeconv"), params = list(output_file = report.html))
}
