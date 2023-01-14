import metacells as mc
import anndata as ad
import numpy as np
import pandas as pd
import seaborn as sb # viz
import matplotlib.pyplot as plt # viz

#  based on the vignette!!!

def clean_genes_and_cells(anndata, properly_sampled_min_cell_total = 800, 
  properly_sampled_max_cell_total = 8000, 
  properly_sampled_max_excluded_genes_fraction = 0.1, 
  exclude_genes = "", exclude_gene_patterns = "", seed = 123456): 
  name = "AnnData"
  
  mc.ut.set_name(anndata, name)
  
  # recommendations from the metacell paper
  excluded_gene_names = list(set(list(exclude_genes) + ['IGHMBP2', 'IGLL1', 'IGLL5', 'IGLON5', 'NEAT1', 'TMSB10', 'TMSB4X']))
  excluded_gene_patterns = list(set(list(exclude_gene_patterns) + ['MT-.*']))
  
  # clean Genes
  mc.pl.analyze_clean_genes(anndata,
                            excluded_gene_names=excluded_gene_names,
                            excluded_gene_patterns=excluded_gene_patterns,
                            random_seed=seed)
  # combine metrics                        
  mc.pl.pick_clean_genes(anndata)
  
  # cleaning cells
  mc.pl.analyze_clean_cells(
      anndata,
      properly_sampled_min_cell_total=properly_sampled_min_cell_total,
      properly_sampled_max_cell_total=properly_sampled_max_cell_total,
      properly_sampled_max_excluded_genes_fraction=properly_sampled_max_excluded_genes_fraction)
  
  mc.pl.pick_clean_cells(anndata)
  
  # extract cleaned data, remove unclean genes and cells
  clean = mc.pl.extract_clean_data(anndata)
  
  return (clean)

def compute_forbidden_genes(clean, suspect_gene_names="", suspect_gene_patterns="", seed = 123456):

  # combine gene names and pattens with recommendations from the metacell vignette                      
  suspect_gene_names = list(set(list(suspect_gene_names) + ['PCNA', 'MKI67', 'TOP2A', 'HIST1H1D',
                        'FOS', 'JUN', 'HSP90AB1', 'HSPA1A', 'ISG15', 'WARS' ]))
  suspect_gene_patterns = list(set(list(suspect_gene_patterns) + [ 'MCM[0-9]', 'SMC[0-9]', 'IFI.*' ]))


  suspect_genes_mask = mc.tl.find_named_genes(clean, names=suspect_gene_names,
                                            patterns=suspect_gene_patterns)

  suspect_gene_names = sorted(clean.var_names[suspect_genes_mask])
  
  mc.pl.relate_genes(clean, random_seed=seed)
  
  module_of_genes = clean.var['related_genes_module']
  suspect_gene_modules = np.unique(module_of_genes[suspect_genes_mask])
  suspect_gene_modules = suspect_gene_modules[suspect_gene_modules >= 0]
  print(suspect_gene_modules)
  
  print (suspect_gene_names)
  
  return (suspect_gene_names)
  
  
  # similarity_of_genes = mc.ut.get_vv_frame(clean, 'related_genes_similarity')
  # for gene_module in suspect_gene_modules:
  #     module_genes_mask = module_of_genes == gene_module
  #     similarity_of_module = similarity_of_genes.loc[module_genes_mask, module_genes_mask]
  #     similarity_of_module.index = \
  #     similarity_of_module.columns = [
  #         '(*) ' + name if name in suspect_gene_names else name
  #         for name in similarity_of_module.index
  #     ]
  #     ax = plt.axes()
  #     sb.heatmap(similarity_of_module, vmin=0, vmax=1, xticklabels=True, yticklabels=True, ax=ax, cmap="YlGnBu")
  #     ax.set_title(f'Gene Module {gene_module}')
  #     plt.show()
  
  
def extract_forbidden_from_modules(clean, forbidden_modules):
  
  # create empty mask 
  forbidden_genes_mask = pd.Series(data = False, index= clean.var_names)
  
  for gene_module in forbidden_modules:
    module_genes_mask = module_of_genes == gene_module
    forbidden_genes_mask |= module_genes_mask
  forbidden_gene_names = sorted(clean.var_names[forbidden_genes_mask])
  print(len(forbidden_gene_names))
  print(' '.join(forbidden_gene_names))
  
  return (clean)
  
def compute_metacells(clean, forbidden_gene_names, seed=123456):
  
  mc.pl.set_max_parallel_piles(mc.pl.guess_max_parallel_piles(clean))
  
  name = "AnnData"
  
  print ("computing metacells...")
  
  # with mc.ut.progress_bar(): # does not work somehow
  mc.pl.divide_and_conquer_pipeline(clean,
                                    forbidden_gene_names=forbidden_gene_names,
                                    #target_metacell_size=...,
                                    random_seed=seed)
                                      
  print ("finished divide and conquer")
                                      
  metacells = mc.pl.collect_metacells(clean, name = name + ".metacells")

  return [clean, metacells]
