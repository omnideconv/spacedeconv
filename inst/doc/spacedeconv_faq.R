## ----eval = FALSE-------------------------------------------------------------
#  rownames(spe) <- rowData(spe)$symbol

## ----eval=FALSE---------------------------------------------------------------
#  spe <- spe[, colData(spe)$in_tissue == TRUE]

## ----availableMethods---------------------------------------------------------
spacedeconv::deconvolution_methods
