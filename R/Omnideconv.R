#' Deconvolute Omnideconv
#' @param spe SpatialExperiment
#' @param signature A signature matrix
#' @param method omnideconv deconvolution method, see omnideconv::deconvolution_methods() for details
#' @param single_cell_object Needed for MOMF and Bisque, either a matrix with cell type annotations in the cell_type_annotations parameter or a singleCellExperiment / AnnData where the cell_type_column_name specifies the corresponding column where the cell type annotation is located in the object
#' @param cell_type_col column name of single_cell_object where cell type information is found
#' @param batch_ids A vector of the ids of samples or individuals
#' @param verbose Whether to produce output on the console
#' @param ... Additional parameters, passed to the selected method

# note, if we decide to treat everything as an object internally we can remove celL_type_annotations from here
deconvolute_omnideconv <- function(spe, signature, method = NULL, single_cell_object = NULL, cell_type_col = "cell_ontology_class", batch_ids = NULL, verbose = FALSE, ...) {
  if (is.null(spe)) {
    stop("Parameter 'spe' is missing or null, but is required.")
  }

  if (is.null(signature)) {
    stop("Parameter 'signature' is missing or null, but is required.")
  }

  if (is.null(method)) {
    stop("Parameter 'method' is missing or null, but is required.")
  }

  if (!(method %in% omnideconv::deconvolution_methods)) {
    stop("Method not supported by omnideconv")
  }

  # extract matrices from object
  bulk_gene_expression <- SingleCellExperiment::counts(spe)

  # extract batch ids???

  # deconvolute
  omnideconv::deconvolute(
    bulk_gene_expression = bulk_gene_expression,
    signature = signature,
    method = method,
    single_cell_object = single_cell_object,
    cell_type_annotations = cell_type_annotations,
    batch_ids = batch_ids,
    cell_type_column_name = cell_type_column_name,
    verbose = verbose,
    ...
  )
}
