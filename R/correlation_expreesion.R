

# build signature matrix from single cell data
signature <- spacedeconv::build_model(
single_cell_obj = single_cell_data_3,
cell_type_col = "celltype_major",
method = "spatialdwls"
)



corr_expr <- function(sig, log = FALSE, cor_method = "pearson"){
  if(log == TRUE){
    sig <- log(counts(sig) +1)
  }

}


