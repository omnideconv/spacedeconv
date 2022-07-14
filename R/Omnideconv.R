#' Build Reference
#'
#' @param single_cell_object SingleCellExperiment
#' @param cell_type_col column containing cell type annotation
#' @param method select deconvolution method
#' @param spatial_object SpatialExperiment
#' @param batch_id_col column of single_cell_object containing batch_ids
#' @param markers List of marker genes, only used for bseqsc
#' @param assay_sc single cell assay to use
#' @param assay_sp spatial assay to use
#' @param ... additional parameters passed to the selected algorithm
#'
#' @returns Gene Expression signature or NULL, if methods don't export signatures
build_model_omnideconv <- function(single_cell_object, cell_type_col = "cell_ontology_class", method = NULL, spatial_object = NULL, batch_id_col = NULL, markers = NULL, assay_sc = "counts", assay_sp= "counts", ...) {
  # TODO checks


  # extract matrices from object
  bulk_gene_expression <- SummarizedExperiment::assay(spatial_object, assay_sp)

  # extract batch ids from single cell object
  batch_ids <- NULL
  if (!is.null(single_cell_object) && !is.null(batch_id_col)) {
    batch_ids <- colData(single_cell_object)[[batch_id_col]]
  }

  signature <- omnideconv::build_model(
    single_cell_object,
    cell_type_column_name = cell_type_col,
    method = method,
    markers = markers,
    bulk_gene_expression = bulk_gene_expression,
    batch_ids = batch_ids,
    ...
  )

  return(signature)
}

#' Deconvolute Omnideconv
#' @param spe SpatialExperiment
#' @param signature A signature matrix
#' @param method omnideconv deconvolution method, see omnideconv::deconvolution_methods() for details
#' @param single_cell_object Needed for MOMF and Bisque, either a matrix with cell type annotations in the cell_type_annotations parameter or a singleCellExperiment / AnnData where the cell_type_column_name specifies the corresponding column where the cell type annotation is located in the object
#' @param cell_type_col column name of single_cell_object where cell type information is found
#' @param batch_id_col column name of single_cell_object where batch ids can be found
#' @param verbose Whether to produce output on the console
#' @param ... Additional parameters, passed to the selected method
deconvolute_omnideconv <- function(spe, signature = NULL, method = NULL, single_cell_object = NULL, cell_type_col = "cell_ontology_class", batch_id_col = NULL, verbose = FALSE, ...) {
  if (is.null(spe)) {
    stop("Parameter 'spe' is missing or null, but is required.")
  }

  if (is.null(method)) {
    stop("Parameter 'method' is missing or null, but is required.")
  }

  if (!(method %in% omnideconv::deconvolution_methods)) {
    stop("Method not supported by omnideconv")
  }

  # extract matrices from object
  bulk_gene_expression <- SingleCellExperiment::counts(spe)

  # extract batch ids from single cell object
  batch_ids <- NULL
  if (!is.null(single_cell_object) && !is.null(batch_id_col)) {
    batch_ids <- colData(single_cell_object)[[batch_id_col]]
  }

  # deconvolute
  omnideconv::deconvolute(
    bulk_gene_expression = bulk_gene_expression,
    signature = signature,
    method = method,
    single_cell_object = single_cell_object,
    batch_ids = batch_ids,
    cell_type_column_name = cell_type_col,
    verbose = verbose,
    ...
  )
}
