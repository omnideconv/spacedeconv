## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----deconvolute, error=FALSE, message=FALSE, warning=FALSE-------------------
library(spacedeconv)
library(SpatialExperiment)

data("single_cell_data_3")
data("spatial_data_3")

## ----normalization, warning=FALSE, eval=FALSE---------------------------------
#  single_cell_data_3 <- spacedeconv::preprocess(single_cell_data_3)
#  spatial_data_3 <- spacedeconv::preprocess(spatial_data_3)
#  
#  single_cell_data_3 <- spacedeconv::normalize(single_cell_data_3, method = "cpm")
#  spatial_data_3 <- spacedeconv::normalize(spatial_data_3, method = "cpm")

## ----build_model, message=FALSE, warning=FALSE, eval=FALSE--------------------
#  signature <- spacedeconv::build_model(
#    single_cell_obj = single_cell_data_3,
#    cell_type_col = "celltype_major",
#    method = "dwls", verbose = T, dwls_method = "mast_optimized", ncores = 10
#  )

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

## -----------------------------------------------------------------------------
spacedeconv::cell_pair_localization(deconv, method = "dwls", cell_type_A = "dwls_T.cells", cell_type_B = "dwls_B.cells", density = TRUE, distance = 0)

## -----------------------------------------------------------------------------
spacedeconv::ripleys_k(deconv, method = "dwls", cell_type = "dwls_B.cells")

