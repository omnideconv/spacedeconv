## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("pak")
#  
#  # minimal installation
#  pak::pkg_install("omnideconv/spacedeconv")
#  
#  # complete installation
#  pak::pkg_install("omnideconv/spacedeconv", dependencies = TRUE)

## ----readData, eval=FALSE-----------------------------------------------------
#  spe <- SpatialExperiment::read10xVisium("path_to_directory")

## ----normalization, eval=FALSE------------------------------------------------
#  spe <- spacedeconv::normalize(spe, method = "cpm")
#  
#  # specify expression assay to use
#  signature <- spacedeconv::build_model(spe,
#    method = "quantiseq",
#    assay_sp = "cpm"
#  )

## ----buildModel, eval=FALSE---------------------------------------------------
#  signature <- spacedeconv::build_model(
#    single_cell_object,
#    cell_type_col = "celltype_major",
#    method = "spotlight",
#    assay_sc = "cpm"
#  )

## ----deconvolution, eval=FALSE------------------------------------------------
#  # save the results to an annotated SpatialExperiment
#  result <- spacedeconv::deconvolute(
#    spatial_object,
#    signature,
#    method = "spotlight"
#  )
#  
#  # return deconvolution results in table form
#  result <- spacedeconv::deconvolute(
#    spatial_object,
#    signature,
#    method = "spotlight",
#    return_object = FALSE
#  )

## ----visualization, eval=FALSE------------------------------------------------
#  # sample does refer to the first column of ColData(spe)
#  # for cell_type input a celltype present in the deconvolution result
#  spacedeconv::plot_celltype(spe, cell_type = "spotlight_B.cells")
#  
#  # umi count
#  spacedeconv::plot_umi_count(spe)

