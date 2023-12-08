## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----deconvolute, error=FALSE, message=FALSE, warning=FALSE-------------------
library(spacedeconv)

data("spatial_data_3")

spe <- preprocess(spatial_data_3)

spe <- spacedeconv::normalize(spe, method = "cpm")

deconv <- spacedeconv::deconvolute(
  spe,
  method = "epic",
  assay_sp = "cpm"
)

## -----------------------------------------------------------------------------
spacedeconv::cell_pair_localization(deconv, method = "epic", cell_type_A = "epic_NK.cell", cell_type_B = "epic_B.cell", density = TRUE, distance = 0)

## -----------------------------------------------------------------------------
spacedeconv::ripleys_k(deconv, method = "epic", cell_type = "epic_B.cell")
