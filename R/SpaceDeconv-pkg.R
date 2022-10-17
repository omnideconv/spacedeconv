#' SpaceDeconv
#'
#' SpaceDeconv Unified R Interface to spatial transcriptomics deconvolution methods
#'
#' @importFrom methods as is
#' @importFrom SingleCellExperiment counts colData
#' @importFrom utils askYesNo install.packages
#' @importFrom remotes install_github
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom ggplot2 aes aes_string discrete_scale element_text labs aes_ ggplot geom_density theme_classic geom_vline scale_fill_viridis_c scale_y_discrete theme annotation_raster geom_sf coord_sf element_blank
#' @importFrom ggridges geom_density_ridges_gradient
#' @importFrom magrittr %>%
#' @importFrom testit assert
#' @importFrom utils read.table
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom S4Vectors DataFrame
#' @importFrom grDevices as.raster
#' @importFrom Matrix colSums
#' @importFrom colorspace scale_fill_continuous_sequential sequential_hcl
#'
#' @name SpaceDeconv-pkg
#' @docType package
NULL
