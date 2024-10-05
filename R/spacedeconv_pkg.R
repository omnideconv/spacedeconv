#' spacedeconv
#'
#' spacedeconv Unified R Interface to spatial transcriptomics deconvolution methods
#'
#' @importFrom methods as is
#' @importFrom SingleCellExperiment counts colData
#' @importFrom SpatialExperiment colData<- spatialCoords
#' @importFrom utils askYesNo install.packages
#' @importFrom remotes install_github
# #' @importFrom Seurat as.SingleCellExperiment
#' @importFrom SummarizedExperiment assay assayNames colData assay<-
#' @importFrom ggplot2 aes aes_string unit element_rect ggtitle coord_fixed geom_text geom_point geom_abline xlab ylab discrete_scale element_text labs aes_ ggplot geom_density theme_classic geom_vline scale_fill_viridis_c scale_y_discrete theme annotation_raster geom_sf coord_sf element_blank scale_fill_brewer scale_fill_manual
#' @importFrom ggridges geom_density_ridges geom_density_ridges_gradient
#' @importFrom magrittr %>%
#' @importFrom testit assert
#' @importFrom ggpubr ggarrange
#' @importFrom utils read.table
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom S4Vectors DataFrame
#' @importFrom grDevices as.raster
#' @importFrom Matrix colSums
#' @importFrom colorspace scale_fill_continuous_sequential sequential_hcl
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom DelayedArray rowSums colSums
#' @importFrom circlize colorRamp2
# #' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom grid gpar
#' @importFrom stats median quantile cor density
#' @importFrom graphics abline par
#' @importFrom grDevices dev.off png
#' @importFrom corrplot corrplot cor.mtest
#' @importFrom multimode locmodes
#' @importFrom methods hasArg is
#' @importFrom decoupleR run_wsum run_aucell run_fgsea run_gsva run_mdt run_mlm run_ora run_udt run_ulm run_viper run_wmean
#' @importFrom psych corr.test
#' @importFrom OmnipathR import_intercell_network
#' @importFrom scales label_number
#' @importFrom Giotto normalizeGiotto
#'
#'
#' @name spacedeconv_pkg
# #' @docType package
NULL
