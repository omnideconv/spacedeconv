# test cell2location

library(reticulate)
reticulate::use_condaenv("cell2loc_env")
reticulate::py_config()
reticulate::py_available()


cl <- import("cell2location")
sc <- import("scanpy")
pd <- import("pandas")
np <- import("numpy")
#ad <- import("anndata") #?

library(anndata) # muss unbedingt geladen werden?


adata_vis = sc$datasets$visium_sge(sample_id="V1_Human_Lymph_Node")


adata_ref = sc$read(
  './data/sc.h5ad',
  backup_url='https://cell2location.cog.sanger.ac.uk/paper/integrated_lymphoid_organ_scrna/RegressionNBV4Torch_57covariates_73260cells_10237genes/sc.h5ad'
)

adata_vis$obs["sample"] <- list(names(adata_vis$uns[["spatial"]])) # add sample information to spots, in case there are more?
adata_vis$var["SYMBOL"] <- adata_vis$var_names # add symbols to "colData"
# das hier oben produziert eine Warning, funktioniert aber trotzdem


# this is the set_index part!!!
# rename genes to ensebml, this could also be symbols
rownames(adata_vis$var) <- adata_vis$var$gene_ids
rownames(adata_ref$var) <- adata_ref$var[["GeneID-2"]]



# this can be extended
selected <- cl$utils$filter_genes(adata_ref)

dev.off() # remove plot?

# subset, unsure here





# MODEL
cl$models$RegressionModel$setup_anndata(adata = adata_ref,
                                        batch_key = "Sample",
                                        labels_key = "Subset",
                                        categorical_covariate_keys = list("Method"))
mod <- cl$models$RegressionModel(adata = adata_ref)

mod$train(max_epochs = 20, use_gpu = FALSE) # use_gpu

mod$plot_history(10)


# export eximated cel abundance = signature
adata_ref_after = mod$export_posterior(adata_ref) # some args missing?


# TODO one can save the model in between
#


inf_aver = adata_ref_after$varm[["means_per_cluster_mu_fg"]]
colnames(inf_aver) <- adata_ref_after$uns[["mod"]][["factor_names"]]  # available cell types



######### now the mapping part

intersect =  np$intersect1d(adata_vis$var_names, rownames(inf_aver) )
intersect = as.vector(intersect)
vin = adata_vis$var_names %in% intersect
# subset both objects

# adata
inf_aver = inf_aver[intersect, ]
adata_vis = adata_vis[ , vin]

cl$models$Cell2location$setup_anndata(adata = adata_vis, batch_key = "sample")

mod = cl$models$Cell2location(adata_vis, cell_state_df = inf_aver, N_cells_per_location = 10, detection_alpha = 200)

