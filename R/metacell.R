library(reticulate)
library(SingleCellExperiment)

mc = reticulate::import("metacells")
ad = reticulate::import("anndata")

seurat = readRDS("~/subsetData/breast_cancer_seurat_object.rds")

sce = Seurat::as.SingleCellExperiment(seurat)


# var = cell  = row
# obs = gene = column
anndata = ad$AnnData(X = SingleCellExperiment::counts(sce),
                     obs = as.data.frame(rowData(sce)), var = as.data.frame(colData(sce)))



# variables
excluded_gene_names = list('IGHMBP2', 'IGLL1', 'IGLL5', 'IGLON5', 'NEAT1', 'TMSB10', 'TMSB4X')
excluded_gene_patterns = list('MT-.*')



mc$pl$analyze_clean_genes(anndata, excluded_gene_names = excluded_gene_names,
                          excluded_gene_patterns = excluded_gene_patterns,
                          random_seed = 123456)

anndata

