## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=FALSE---------------------------------------------------------------
library(spacedeconv)
library(SpatialExperiment)
spacedeconv::deconvolution_methods

## ----loadData, message=FALSE, eval = TRUE-------------------------------------
data("single_cell_data_3")
data("spatial_data_3")

## ----viewData, eval=FALSE-----------------------------------------------------
#  single_cell_data_3
#
#  # for performance reasons we are subsampling the data
#
#  # single_cell_data_3 <- subsetSCE(single_cell_data_3, cell_type_col = "celltype_major", ncells = 180)

## ----normalization, warning=FALSE, eval=FALSE---------------------------------
#  single_cell_data_3 <- spacedeconv::preprocess(single_cell_data_3)
#  spatial_data_3 <- spacedeconv::preprocess(spatial_data_3)
#
#  single_cell_data_3 <- spacedeconv::normalize(single_cell_data_3, method = "cpm")
#  spatial_data_3 <- spacedeconv::normalize(spatial_data_3, method = "cpm")

## ----viewData3----------------------------------------------------------------
names(colData(single_cell_data_3))

## ----build_model, message=FALSE, warning=FALSE, eval=FALSE--------------------
#  signature <- spacedeconv::build_model(
#    single_cell_obj = single_cell_data_3,
#    cell_type_col = "celltype_major",
#    method = "dwls", verbose = T, dwls_method = "mast_optimized", ncores = 10
#  )

## -----------------------------------------------------------------------------
signature <- readRDS(system.file("extdata", "signature_dwls.rds", package = "spacedeconv"))

## ----eval=TRUE----------------------------------------------------------------
knitr::kable(round(signature[1:10, ], 4))

## ----deconvolution, eval=FALSE, message=FALSE, warning=FALSE------------------
#  deconv <- spacedeconv::deconvolute(
#    spatial_obj = spatial_data_3,
#    single_cell_obj = single_cell_data_3,
#    cell_type_col = "celltype_major",
#    method = "dwls",
#    signature = signature,
#    assay_sp = "cpm"
#  )

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
deconv <- readRDS(system.file("extdata", "deconv_dwls.rds", package = "spacedeconv"))

## ----accessColData------------------------------------------------------------
available_results(deconv)

## ----visualization, fig.width=12, fig.height=12-------------------------------
# plot all available results: provide the tool as parameter
spacedeconv::plot_celltype(deconv,
  cell_type = "dwls_B.cells",
  density = FALSE,
  smooth = T,
  title_size = 12
)

## ----visualization2, fig.width=12, fig.height=12------------------------------
# ... or plot a specific result
# spacedeconv::plot_celltype(deconv,
#   cell_type = "card_Cancer.Epithelial",
#   density = FALSE,
#   smooth = T,
#   title_size = 12
# )
