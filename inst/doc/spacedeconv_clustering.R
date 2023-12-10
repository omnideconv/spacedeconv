## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(spacedeconv)
library(SpatialExperiment)

data("single_cell_data_3")
data("spatial_data_3")

## ----deconv, eval=FALSE-------------------------------------------------------
#  single_cell_data_3 <- spacedeconv::preprocess(single_cell_data_3)
#  spatial_data_3 <- spacedeconv::preprocess(spatial_data_3)
#  
#  single_cell_data_3 <- spacedeconv::normalize(single_cell_data_3, method = "cpm")
#  spatial_data_3 <- spacedeconv::normalize(spatial_data_3, method = "cpm")
#  deconv <- deconvolute(spe, method = "epic", assay_sc = "cpm")

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

## ----eval=FALSE---------------------------------------------------------------
#  cluster <- spacedeconv::cluster(deconv, data = "expression", clusres = 0.5)

## -----------------------------------------------------------------------------
cluster <- readRDS(system.file("extdata", "cluster.rds", package = "spacedeconv"))

## ----clusterExpression, eval=FALSE--------------------------------------------
#  plot_celltype(cluster, "cluster", density = F) # plot the clustering stored in this object

## ----topFeatures--------------------------------------------------------------
get_cluster_features(cluster, clusterid = "cluster_expression_0.5", spmethod = "dwls")

