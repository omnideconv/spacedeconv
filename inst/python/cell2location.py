import scanpy as sc
import anndata
import pandas as pd
import cell2location
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#from matplotlib import rcParams
#rcParams["pdf.fonttype"] = 42


#sp_obj = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
#sp_obj.obs['sample'] = list(sp_obj.uns['spatial'].keys())[0] # add sample information

#sp_obj.var['SYMBOL'] = sp_obj.var_names # add rownames as column, in this case
#sp_obj.var.set_index('gene_ids', drop=True, inplace=True) # ensembl as rownames???

# adata_ref = sc.read(
#     f'./data/sc.h5ad',
#    backup_url='https://cell2location.cog.sanger.ac.uk/paper/integrated_lymphoid_organ_scrna/RegressionNBV4Torch_57covariates_73260cells_10237genes/sc.h5ad'
# )
# 
# adata_ref.var['SYMBOL'] = adata_ref.var.index
# rename 'GeneID-2' as necessary for your data
# adata_ref.var.set_index('GeneID-2', drop=True, inplace=True)
# 
# delete unnecessary raw slot (to be removed in a future version of the tutorial)
# del adata_ref.raw

def py_build_model_cell2location(adata_ref, 
                                  epochs = 20, 
                                  sample = "Sample", 
                                  cell_type_column = "celltype_major", 
                                  cell_count_cutoff=5, 
                                  cell_percentage_cutoff=0.03, 
                                  nonz_mean_cutoff=1.12, gpu = False):
                                  
  """
  Build a model using cell2location
  
  Parameters
  ----------
  sc_object: anndata
    AnnData Object containing single-cell reference gene expression 
  sp_object: anndata
    AnnData Object containing spatial expression data
  epochs: integer
    Number of epochs to train the model
  cell_count_cutoff: integer
    Cell2location parameter
  cell_percentage_cutoff: float
    Cell2location parameter
  nonz_mean_cutoff: float
    cell2location parameter
  """
  
  print ("THis is the type", type(adata_ref))
  
  # a filtering step, outsource parameters
  from cell2location.utils.filtering import filter_genes
  # filter genes#
  mpl.pyplot.ioff()
  selected = filter_genes(adata_ref, 
    cell_count_cutoff = cell_count_cutoff, 
    cell_percentage_cutoff2 = cell_percentage_cutoff, 
    nonz_mean_cutoff = nonz_mean_cutoff)
  adata_ref = adata_ref[:, selected].copy() # filter

  cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                          # 10X reaction / sample / batch
                          batch_key=sample,
                          # cell type, covariate used for constructing signatures
                          labels_key=cell_type_column
                         )
                        
  from cell2location.models import RegressionModel
  mod = RegressionModel(adata_ref)
  print ("Training Model")
  mod.train(max_epochs = epochs)
  
  # export model to anndata object
  print ("finished training, extracting results")
  adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False}
  )

  # export estimated expression in each cluster
  if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
      signature = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                      for i in adata_ref.uns['mod']['factor_names']]].copy()
  else:
      signature = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                      for i in adata_ref.uns['mod']['factor_names']]].copy()
                                      
  # rename celltypes
  signature.columns = adata_ref.uns['mod']['factor_names']

  return signature

#################
# DECONVOLUTION #
#################
def py_deconvolute_cell2location(sp_obj, signature, epochs= 30000, n_cell=10, alpha=20, gpu = True):
  """
  Deconvolute using cell2location 
  
  Parameters
  ----------
  sp_obj: AnnData
    AnnData containing spatial gene expression 
  signature: DataFrame
    DataFrame with gene expression signatures
  epochs: integer
    number of epochs used for trainging
  n_cell: Integer
    Number of expected cells per spot
  alpha: float
    detection alpha
  gpu: boolean
    use gpu acceleration
  """

  # spatial mapping 
  # find shared genes and subset both anndata and reference signatures
  intersect = np.intersect1d(sp_obj.var_names, signature.index)
  sp_obj = sp_obj[:, intersect].copy()
  signature = signature.loc[intersect, :].copy()
  
  # prepare anndata for cell2location model
  cell2location.models.Cell2location.setup_anndata(adata=sp_obj, batch_key="sample_id") ############ PARAMETER
  
  # create and train the model
  mod = cell2location.models.Cell2location(
      sp_obj, cell_state_df=signature,
      # the expected average cell abundance: tissue-dependent
      # hyper-prior which can be estimated from paired histology:
      N_cells_per_location=n_cell,
      # hyperparameter controlling normalisation of
      # within-experiment variation in RNA detection:
      detection_alpha=alpha
  )
  
  # train model 
  mod.train(max_epochs=100,
            # train using full data (batch_size=None)
            batch_size=None,
            # use all data points in training because
            # we need to estimate cell abundance at all locations
            train_size=1,
            use_gpu=gpu)
            
  # export fractions
  sp_obj = mod.export_posterior(
    sp_obj, sample_kwargs = {'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': gpu}
  )
  
  # rename
  sp_obj.obs[sp_obj.uns['mod']['factor_names']] = sp_obj.obsm['q05_cell_abundance_w_sf']
  
  return sp_obj.obs[sp_obj.uns["mod"]["factor_names"]]
