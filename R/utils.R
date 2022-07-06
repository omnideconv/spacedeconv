#' get matrices from SingleCellExperiment
#' @param single_cell_object SingleCellExperiment
#' @param cell_type_col column containing the cell type
#' @returns list containing the expression matrix and cell type annotation vector

getMatricesFromSCE <- function(single_cell_object, cell_type_col = "cell_ontology_class") {
  # tests

  # test if object valid

  # test if cell_type_col in names(colData())

  counts <- as(counts(single_cell_object), "dgCMatrix") # count matrix as sparse matrix
  cell_type_annotation <- as.character(colData(single_cell_object)[[cell_type_col]])


  return(list(counts = counts, cell_type_annotation = cell_type_annotation))
}

#' Check if Column exists in object

checkCol <- function(object, column){
  return (column %in% names(colData(object)))
}
