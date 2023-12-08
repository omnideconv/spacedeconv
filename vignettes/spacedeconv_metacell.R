## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(spacedeconv)
library(SpatialExperiment)
data("single_cell_data_1")
single_cell_data_1
ad <- spe_to_ad(single_cell_data_1) # convert to anndata

## ----clean--------------------------------------------------------------------
filtered <- clean_genes_and_cells(ad)

## -----------------------------------------------------------------------------
suspect_genes <- compute_forbidden_genes(filtered)

## -----------------------------------------------------------------------------
metacells <- compute_metacells(filtered, suspect_genes,
  cell_type_col = "celltype_major",
  abundance_score = 0.9
)

## -----------------------------------------------------------------------------
metacells
