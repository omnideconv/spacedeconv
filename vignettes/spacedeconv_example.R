## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=FALSE---------------------------------------------------------------
library(spacedeconv)
spacedeconv::deconvolution_methods

## ----loadData, message=FALSE--------------------------------------------------
library(spacedeconv)
library(SpatialExperiment)

data("single_cell_data_3")
data("spatial_data_3")

## ----viewData-----------------------------------------------------------------
single_cell_data_3

# for performance reasons we are subsampling the data

single_cell_data_3 <- subsetSCE(single_cell_data_3, cell_type_col = "celltype_major", ncells = 180)

## ----normalization, warning=FALSE---------------------------------------------
single_cell_data_3 <- spacedeconv::preprocess(single_cell_data_3)
spatial_data_3 <- spacedeconv::preprocess(spatial_data_3)

single_cell_data_3 <- spacedeconv::normalize(single_cell_data_3, method = "cpm")
spatial_data_3 <- spacedeconv::normalize(spatial_data_3, method = "cpm")

## ----viewData3----------------------------------------------------------------
names(colData(single_cell_data_3))

## ----build_model, message=FALSE, warning=FALSE, eval=FALSE--------------------
#  signature <- spacedeconv::build_model(
#    single_cell_obj = single_cell_data_3,
#    cell_type_col = "celltype_major",
#    method = "card"
#  )

## ----eval=FALSE---------------------------------------------------------------
#  knitr::kable(round(signature[1:10, ], 4))

## ----deconvolution, message=FALSE, warning=FALSE------------------------------
deconv <- spacedeconv::deconvolute(
  spatial_obj = spatial_data_3,
  single_cell_obj = single_cell_data_3,
  cell_type_col = "celltype_major",
  method = "card",
  batch_id_col = "orig.ident",
  assay_sp = "cpm",
  assay_sc = "cpm"
)

## ----accessColData------------------------------------------------------------
available_results(deconv)

