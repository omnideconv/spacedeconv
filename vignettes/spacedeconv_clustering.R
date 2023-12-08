## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(spacedeconv)
data("spatial_data_3")

## ----deconv-------------------------------------------------------------------
spe <- preprocess(spatial_data_3)
spe <- normalize(spe, method = "cpm")
deconv <- deconvolute(spe, method = "epic", assay_sc = "cpm")

