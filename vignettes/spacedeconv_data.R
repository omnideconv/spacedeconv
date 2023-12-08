## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----loadingData, eval=FALSE--------------------------------------------------
#  library(spacedeconv)
#  data("single_cell_data_1")
#  data("single_cell_data_2")
#  data("single_cell_data_3")
#  data("single_cell_data_4")
#  data("spatial_data_1")
#  data("spatial_data_2")
#  data("spatial_data_3")
#  data("spatial_data_4")

## ----annotationParameters, eval=FALSE-----------------------------------------
#  signature <- build_model(single_cell_data_1,
#    method = "dwls",
#    cell_type_col = "celltype_major",
#  )
#  
#  # some methods require batch_id information as well
#  sigature <- build_model(single_cell_data_1,
#    method = "scdc",
#    cell_type_col = "celltype_major",
#    batch_id_col = "orig.ident"
#  )

## ----loadSpatialExperiment, eval=FALSE----------------------------------------
#  spe <- SpatialExperiment::read10xVisium(path_to_directory)

## ----seuratConversion, eval=FALSE---------------------------------------------
#  spe <- seurat_to_spatialexperiment(seurat_object)

## ----deconvolution, eval = FALSE----------------------------------------------
#  spe <- normalize(spe, method = "cpm")
#  
#  # make sure to use cpm assay in deconvolution step
#  deconvolution <- deconvolute(spe, method = "quantiseq", assay_sp = "cpm")

