#' Build Reference
#'
#' @param single_cell_obj SingleCellExperiment
#' @param cell_type_col column containing cell type annotation
#' @param method select deconvolution method
#' @param spatial_obj SpatialExperiment
#' @param batch_id_col column of single_cell_obj containing batch_ids
#' @param markers List of marker genes, only used for bseqsc
#' @param assay_sc single cell assay to use
#' @param assay_sp spatial assay to use
#' @param ... additional parameters passed to the selected algorithm
#'
#' @returns Gene Expression signature or NULL, if methods don't export signatures
build_model_omnideconv <- function(single_cell_obj, cell_type_col = "cell_ontology_class", method = NULL, spatial_obj = NULL, batch_id_col = NULL, markers = NULL, assay_sc = "counts", assay_sp= "counts", ...) {
  # TODO checks

  # check if requested assay exists
  if (!assay_sc %in% names(SummarizedExperiment::assays(single_cell_obj))) {
    message(
      "requested assay ",
      assay_sc,
      " not available in expression object. Using first available assay."
    )
    assay_sc <- names(SummarizedExperiment::assays(single_cell_obj))[1] # change to first available assay request not available
  }

  # check if requested assay exists
  if (!is.null(spatial_obj) && !assay_sp %in% names(SummarizedExperiment::assays(spatial_obj))) {
    message(
      "requested assay ",
      assay_sp,
      " not available in expression object. Using first available assay."
    )
    assay_sp <- names(SummarizedExperiment::assays(spatial_obj))[1] # change to first available assay request not available
  }

  # extract matrices from object
  bulk_gene_expression <- NULL
  if (!is.null(spatial_obj)){
    bulk_gene_expression <- SummarizedExperiment::assay(spatial_obj, assay_sp)
  }

  # extract batch ids from single cell object
  batch_ids <- NULL
  if (!is.null(single_cell_obj) && !is.null(batch_id_col)) {
    batch_ids <- colData(single_cell_obj)[[batch_id_col]]
  }

  signature <- omnideconv::build_model(
    single_cell_obj,
    cell_type_column_name = cell_type_col,
    method = method,
    markers = markers,
    bulk_gene_expression = bulk_gene_expression,
    batch_ids = batch_ids,
    assay_name = assay_sc,
    ...
  )

  return(signature)
}

#' Deconvolute Omnideconv
#' @param spatial_obj SpatialExperiment
#' @param signature A signature matrix
#' @param method omnideconv deconvolution method, see omnideconv::deconvolution_methods() for details
#' @param single_cell_obj Needed for MOMF and Bisque, either a matrix with cell type annotations in the cell_type_annotations parameter or a singleCellExperiment / AnnData where the cell_type_column_name specifies the corresponding column where the cell type annotation is located in the object
#' @param cell_type_col column name of single_cell_obj where cell type information is found
#' @param batch_id_col column name of single_cell_obj where batch ids can be found
#' @param assay_sc single cell assay to use
#' @param assay_sp spatial assay to use
#' @param verbose Whether to produce output on the console
#' @param ... Additional parameters, passed to the selected method
deconvolute_omnideconv <- function(spatial_obj, signature = NULL, method = NULL, single_cell_obj = NULL, cell_type_col = "cell_ontology_class", batch_id_col = NULL, assay_sc = "counts", assay_sp="counts", verbose = FALSE, ...) {
  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is missing or null, but is required.")
  }

  if (is.null(method)) {
    stop("Parameter 'method' is missing or null, but is required.")
  }

  if (!(method %in% omnideconv::deconvolution_methods)) {
    stop("Method not supported by omnideconv")
  }

  # extract matrices from object
  bulk_gene_expression <- NULL
  if (!is.null(spatial_obj)){
    bulk_gene_expression <- SummarizedExperiment::assay(spatial_obj, assay_sp)
  }

  # extract batch ids from single cell object
  batch_ids <- NULL
  if (!is.null(single_cell_obj) && !is.null(batch_id_col)) {
    batch_ids <- colData(single_cell_obj)[[batch_id_col]]
  }

  # deconvolute
  omnideconv::deconvolute(
    bulk_gene_expression = bulk_gene_expression,
    signature = signature,
    method = method,
    single_cell_object = single_cell_obj,
    batch_ids = batch_ids,
    cell_type_column_name = cell_type_col,
    verbose = verbose,
    assay_name = assay_sc,
    ...
  )
}
