#' No signature calculated, just call the deconvolute method
#' @returns NULL

build_model_immunedeconv <- function() {
  message("This method does not build a signature, just call the deconvolute method")
  return(NULL)
}

#' Deconvolute Immundeconv
#' @param spe SpatialExperiment
#' @param method deconvolution algorithm
#' @param assay_sp assay of spatial object to use
#' @param ... further parameters passed to the selected method
deconvolute_immunedeconv <- function(spe, method = NULL, assay_sp = "counts", ...) {
  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if (is.null(method)) {
    stop("Parameter 'method' is missing or null, but is required")
  }

  # manual workaround for xCell
  if (method == "xcell" && requireNamespace("xCell")){
    xCell.data <- xCell::xCell.data
  }

  # check for HGNC Symbols!!!
  rownames(spe) <- toupper(rownames(spe)) # temporary fix

  # check if requested assay exists
  if (!assay_sp %in% names(SummarizedExperiment::assays(spe))) {
    message(
      "requested assay ",
      assay_sp,
      " not available in expression object. Using first available assay."
    )
    assay_sp <- names(SummarizedExperiment::assays(spe))[1] # change to first available assay request not available
  }

  # extract bulk
  bulk <- NULL
  if (!is.null(spe)){
    bulk <- SummarizedExperiment::assay(spe, assay_sp)
    bulk <- as.matrix(bulk)
  }

  # if method is "cibersort" remove spots with zero expression in signature genes
  if (method=="cibersort"){
    # check if cibersort variables are set

    testit::assert("CIBERSORT.R is provided", exists("cibersort_binary", envir = config_env))
    testit::assert("CIBERSORT signature matrix is provided", exists("cibersort_mat", envir = config_env))

    sig = read.table(get("cibersort_mat", envir = config_env), header = TRUE, sep="\t", row.names = 1, check.names = FALSE)

    #intersect genes
    sig_genes <- rownames(sig)
    bulk_genes <- rownames(bulk)
    bulk_in_sig <- bulk_genes %in% sig_genes
    tmp <- colSums(bulk[bulk_in_sig, ])==0
    if (sum(tmp)>0){
      message("Removing unsusable spots")
      bulk = bulk[, !tmp]
    }
  }

  # some parameters are handled with the ... option
  deconv <- immunedeconv::deconvolute(
    gene_expression = bulk,
    method = method,
    indications = NULL,
    tumor = TRUE,
    arrays = FALSE
  )

  return(convertImmunedeconvMatrix(deconv))
}


#' Deconvolute Immunedeconv mouse
#' @param spe SpatialExperiment
#' @param method deconvolution algorithm
#' @param rmgenes genes to remove from the analysis
#' @param algorithm statistical algorithm for SeqImmuCC (ignored by all other methods)
#' @param assay_sp assay of spatial object to use
#' @param ... additional parameters passed to function
#'
deconvolute_immunedeconv_mouse <- function(spe, method = NULL, rmgenes = NULL, algorithm = NULL, assay_sp = "counts", ...) {
  if (is.null(spe)) {
    stop("Parameter 'spe' is null or missing, but is required")
  }

  if (is.null(method)) {
    stop("Parameter 'method' is missing or null, but is required")
  }

  # check if requested assay exists
  if (!assay_sp %in% names(SummarizedExperiment::assays(spe))) {
    message(
      "requested assay ",
      assay_sp,
      " not available in single cell object. Using first available assay."
    )
    assay_sp <- names(SummarizedExperiment::assays(spe))[1] # change to first available assay request not available
  }

  # extract bulk
  bulk <- NULL
  if (!is.null(spe)){
    bulk <- SummarizedExperiment::assay(spe, assay_sp)
    bulk <- as.matrix(bulk)
  }

  deconv <- immunedeconv::deconvolute_mouse(
    gene_expression_matrix = bulk,
    method = method,
    rmgenes = rmgenes,
    algorithm = algorithm
  )

  return(convertImmunedeconvMatrix(deconv))
}

#' Convert Immunedeconv Matrix to match SpaceDeconv Format
#' @param deconvResult immunedeconv result matrix
#' @returns transformed matrix where cell types are columns
convertImmunedeconvMatrix <- function(deconvResult) {
  # TODO Checks


  # use column 1 as rownames and transpose
  result <- deconvResult[, 2:ncol(deconvResult)]
  result <- as.matrix(result)
  rownames(result) <- deconvResult$cell_type
  result <- t(result)

  return(result)
}
