#' Plot the number of genes per metacells
#'
#' @param metacell SingleCellExperiment containing metacells
#'
#' @export
plot_metacell_nCells <- function(metacell) {
  plot(density(colData(metacell)$grouped), main = "Number of Cells per metacell")
}

#' Plot metacell abundance score
#'
#' @param metacell SingleCellExperiment containing metacells
#'
#' @export
plot_metacell_abundance <- function(metacell) {
  plot(density(colData(metacell)$percentage), main = "Abundance Score per metacell")
}
