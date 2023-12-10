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

## ----clean, eval=FALSE--------------------------------------------------------
#  filtered <- clean_genes_and_cells(ad)

## ----eval=FALSE---------------------------------------------------------------
#  suspect_genes <- compute_forbidden_genes(filtered)

## ----eval=FALSE---------------------------------------------------------------
#  metacells <- compute_metacells(filtered, suspect_genes,
#    cell_type_col = "celltype_major",
#    abundance_score = 0.9
#  )

## -----------------------------------------------------------------------------
metacells <- readRDS(system.file("extdata", "metacells.rds", package = "spacedeconv"))

## -----------------------------------------------------------------------------
metacells

