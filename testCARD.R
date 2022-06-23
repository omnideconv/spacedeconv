# test CARD

library(CARD)
library(SpatialExperiment)
library(magrittr)

sp <- TENxVisiumData::MouseKidneyCoronal()
rownames(sp) <- rowData(sp)$symbol # overwrite EnsemblID

sc <- TabulaMurisSenisData::TabulaMurisSenisDroplet(tissues= "Kidney")$Kidney
sc <- sc[, sc$age=="18m"]
sc <- sc[, !sc$free_annotation %in% c("nan", "CD45")] # remove false data

# downsampling for performance reasons, copied from MarcElosua/SPOTlight
idx <- split(seq(ncol(sc)), sc$free_annotation)
n_cells <- 50
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})
sc <- sc[, unlist(cs_keep)]






spExp <- counts(sp) # expression
spCoord <- SpatialExperiment::spatialCoords(sp) %>% as.data.frame() # coords
colnames(spCoord) <- c("x", "y") # need to be x and y

scExp <- counts(sc) %>% as("dgCMatrix")




scMeta <- data.frame(cellType = colData(sc)$cell_ontology_class, sampleInfo=colData(sc)$mouse.id, row.names = colnames(sc))


cardObj <- createCARDObject(sc_count = scExp, sc_meta = scMeta, spatial_count = spExp, spatial_location = spCoord, sample.varname = "sampleInfo", ct.varname = "cellType", ct.select = cts)

cardObj <- CARD_deconvolution(cardObj)
