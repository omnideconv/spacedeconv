library("Giotto")
library("SingleCellExperiment")


path_file <- "/home/czackl/Dokumente/exchange/BREAST_CANCER_BLOCK_A_SECTION_1"
path_result <- "/home/czackl/Dokumente/exchange/"
# path_file <- "/home/czackl/Documents/DataSpaceDeconv/BREAST_CANCER_BLOCK_A_SECTION_1"
# path_result <- "/home/czackl/Documents/"
# path_python <- reticulate::py_config()$python # python path


# create Giotto Instructions
instrs <- Giotto::createGiottoInstructions(
  save_dir = path_result,
  save_plot = TRUE,
  show_plot = TRUE
)

# create Giotto Object
visium <- Giotto::createGiottoVisiumObject(
  visium_dir = path_file, expr_data = "raw",
  png_name = "tissue_lowres_image.png", instructions = instrs
)


# position image properly
visium <- updateGiottoImage(visium, image_name = "image", xmax_adj = 3500, xmin_adj = 3000, ymax_adj = 4800, ymin_adj = 2000)

# use only in_tissue cells
metadata <- pDataDT(visium)
in_tissue_bc <- metadata[in_tissue == 1]$cell_ID
visium <- subsetGiotto(visium, cell_ids = in_tissue_bc)

# optional filtering?????
visium <- filterGiotto(visium, expression_threshold = 1, verbose = T) # remove all zero genes


visium <- normalizeGiotto(visium)

# create Reference Matrix

sc <- TabulaMurisSenisData::TabulaMurisSenisDroplet(tissues = "Kidney")$Kidney

mat <- counts(sc)
cell_types <- as.vector(sc$cell_ontology_class)
signature <- Giotto::makeSignMatrixDWLSfromMatrix(mat, sign_gene = rownames(sc), cell_type_vector = cell_types)
# signature <- signature[1:100, ] # subset for performance
signature <- signature[rowSums(signature) != 0, ] # remove all zero genes?????


library("org.Hs.eg.db")
ids <- mapIds(org.Hs.eg.db, keys = toupper(rownames(signature)), keytype = "SYMBOL", column = "ENSEMBL") # translate symbol to ensemblID
rownames(signature) <- ids # set rownames
signature <- signature[!is.na(rownames(signature)), ] # remove all unmapped genes, is there a better way to do this?
signature <- signature[rowSums(signature) != 0, ] # remove all zero genes, eigen error!


# head(signature)


# visium <- runHyperGeometricEnrich(visium, signature)


# i dont know why but we need to cluster
visium <- calculateHVG(visium)
visium <- runPCA(visium)
visium <- createNearestNetwork(visium)
visium <- doLeidenCluster(visium, name = "leiden_clus")

message("Started Deconvolution")

deconvolution <- runDWLSDeconv(visium, sign_matrix = signature, n_cell = 10, return_gobject = FALSE)

saveRDS(deconvolution, "/home/czackl/Dokumente/exchange/result2.rds")
