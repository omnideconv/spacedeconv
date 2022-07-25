# test cell2location

library(reticulate)
reticulate::use_condaenv("cell2loc_env")
reticulate::py_config()
reticulate::py_available()



cl <- import("cell2location")
sc <- import("scanpy")
ad <- import("anndata")
pd <- import("pandas")


adata_vis = sc$datasets$visium_sge(sample_id="V1_Human_Lymph_Node")


adata_ref = sc$read(
  './data/sc.h5ad',
  backup_url='https://cell2location.cog.sanger.ac.uk/paper/integrated_lymphoid_organ_scrna/RegressionNBV4Torch_57covariates_73260cells_10237genes/sc.h5ad'
)

adata_vis$obs["sample"] <- list(names(adata_vis$uns[["spatial"]])) # add sample information to spots, in case there are more?
adata_vis$var["SYMBOL"] <- adata_vis$var_names # add symbols to "colData"


# this is the set_index part!!!
# rename genes to ensebml
rownames(adata_vis$var) <- adata_vis$var$gene_ids
rownames(adata_ref$var) <- adata_ref$var[["GeneID-2"]]



# this can be extended
selected <- cl$utils$filter_genes(adata_ref)


# subset, unsure here





# MODEL
cl$models$RegressionModel$setup_anndata(adata = adata_ref,
                                        batch_key = "Sample",
                                        labels_key = "Subset",
                                        categorical_covariate_keys = list("Method"))
mod <- cl$models$RegressionModel(adata = adata_ref)

mod$train(max_epochs = 2) # use_gpu



