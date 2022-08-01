#' List of supported deconvolution methods
#' @export
#'
deconvolution_methods <- c(
  # spatial
  "RCTD" = "rctd",
  "SPOTlight" = "spotlight",
  "CARD" = "card",
  "spatialDWLS" = "spatialdwls",
  "cell2location" = "cell2location",
  "SpatialDecon" = "spatialdecon",
  "DestVi" = "destvi",
  # omnideconv
  "AutoGeneS" = "autogenes",
  "Bisque" = "bisque",
  "BSeq-sc" = "bseqsc",
  "CIBERSORTx" = "cibersortx",
  "CDSeq" = "cdseq",
  "CPM" = "cpm",
  "DWLS" = "dwls",
  "MOMF" = "momf",
  "MuSiC" = "music",
  "Scaden" = "scaden",
  "SCDC" = "scdc",
  # immunedeconv
  "MCPcounter" = "mcp_counter",
  "EPIC" = "epic",
  "quanTIseq" = "quantiseq",
  "xCell" = "xcell",
  "CIBERSORT" = "cibersort",
  "CIBERSORT (abs.)" = "cibersort_abs",
  "TIMER" = "timer",
  "ConsensurTME" = "consensus_tme",
  "ABIS" = "abis",
  "ESTIMATE" = "estimate",
  # immunedeconv mouse
  "mMCPcounter" = "mmcp_counter",
  "seqImmuCC" = "seqimmucc",
  "DCQ" = "dcq",
  "BASE" = "base"
)


#' Build Reference
#' @param single_cell_object Single Cell Expression Object
#' @param cell_type_col Name of the anntotation column containing cell type information
#' @param method signature calculation method
#' @param verbose Display more information on console
#' @param spatial_object SpatialExperiment, required for SPOTlight
#' @param batch_id_col column of singleCellExperiment containing batch ids
#' @param assay_sc assay of single cell object to use
#' @param assay_sp assay of spatial object to use
#' @param ... additional parameters passed to the functions
#' @export
build_model <- function(single_cell_object, cell_type_col = "cell_ontology_class", method = NULL, verbose = FALSE, spatial_object = NULL, batch_id_col = NULL, assay_sc = "counts", assay_sp = "counts", ...) {
  if (is.null(single_cell_object)) {
    stop("Parameter 'single_cell_object' missing or null, but is required")
  }

  # check if method null or not supported
  if (is.null(method) || !(method %in% deconvolution_methods)) {
    stop("Parameter 'method' is null or not supported")
  }

  # if got the methods name and not the token
  if (method %in% names(deconvolution_methods)) {
    method <- deconvolution_methods[[method]]
  }
  method <- tolower(method)
  check_and_install(method)

  # convert data
  single_cell_object <- convert_to_sce(single_cell_object)


  # check if cell_type_col in names(colData(sc)^)


  signature <- switch(method,
    rctd = {
      build_model_rctd()
    },
    spotlight = {
      build_model_spotlight(single_cell_object, cell_type_col, spatial_object, assay_sc = assay_sc, assay_sp = assay_sp, ...)
    },
    card = {
      build_model_card()
    },
    spatialdwls = {build_model_spatial_dwls(sce, assay_sc = assay_sc, marker_method = "scran", cell_type_col = cell_type_col, ...)},

    ##############
    # omnideconv #
    ##############
    autogenes = {
      build_model_omnideconv(single_cell_object, cell_type_col, method = "autogenes", spatial_object, batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose, ...)
    },
    bisque = {
      build_model_omnideconv(single_cell_object, cell_type_col, method = "bisque", spatial_object, batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose, ...)
    },
    bseqsc = {
      build_model_omnideconv(single_cell_object, cell_type_col, method = "bseqsc", spatial_object, batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose, ...)
    },
    cibersortx = {
      build_model_omnideconv(single_cell_object, cell_type_col, method = "cibersortx", spatial_object, batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose, ...)
    },
    cdseq = {
      build_model_omnideconv(single_cell_object, cell_type_col, method = "cdseq", spatial_object, batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose, ...)
    },
    cpm = {
      build_model_omnideconv(single_cell_object, cell_type_col, method = "cpm", spatial_object, batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose, ...)
    },
    dwls = {
      build_model_omnideconv(single_cell_object, cell_type_col, method = "dwls", spatial_object, batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose, ...)
    },
    momf = {
      build_model_omnideconv(single_cell_object, cell_type_col, method = "momf", spatial_object, batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose, ...)
    },
    music = {
      build_model_omnideconv(single_cell_object, cell_type_col, method = "music", spatial_object, batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose, ...)
    },
    scaden = {
      build_model_omnideconv(single_cell_object, cell_type_col, method = "scaden", spatial_object, batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose, ...)
    },
    scdc = {
      build_model_omnideconv(single_cell_object, cell_type_col, method = "scdc", spatial_object, batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose, ...)
    },

    ################
    # immunedeconv #
    ################
    mcp_counter = {
      build_model_immunedeconv()
    },
    epic = {
      build_model_immunedeconv()
    },
    quantiseq = {
      build_model_immunedeconv()
    },
    xcell = {
      build_model_immunedeconv()
    },
    cibersort = {
      build_model_immunedeconv()
    },
    cibersort_abs = {
      build_model_immunedeconv()
    },
    timer = {
      build_model_immunedeconv()
    },
    consensus_tme = {
      build_model_immunedeconv()
    },
    abis = {
      build_model_immunedeconv()
    },
    estimate = {
      build_model_immunedeconv()
    },

    ######################
    # immunedeconv mouse #
    ######################
    mmcp_counter = {
      build_model_immunedeconv()
    },
    seqimmucc = {
      build_model_immunedeconv()
    },
    dcq = {
      build_model_immunedeconv()
    },
    base = {
      build_model_immunedeconv()
    }
  )

  return(signature)
}


#' Deconvolution
#' @param spatial_object A SpatialExperiment
#' @param signature Gene Expression Signature
#' @param single_cell_object A SingleCellExperiment
#' @param cell_type_col Column name of the single_cell_object where the cell type can be found
#' @param method Deconvolution Method to use, see deconvolution_methods() for a full list of available methods
#' @param batch_id_col column name of batch ids in single cell object
#' @param assay_sc which single cell assay to use for computation
#' @param assay_sp which spatial assay to use for computation
#' @param return_object Return an Object or result Table, TRUE = Object
#' @param verbose display more detailed information
#' @param ... Further parameters passed to the selected deconvolution method
#' @returns The deconvolution result as a table
#' @export
deconvolute <- function(spatial_object, signature = NULL, single_cell_object = NULL, cell_type_col = "cell_ontology_class", method = NULL, batch_id_col = NULL, assay_sc = "counts", assay_sp = "counts", return_object = TRUE, verbose = FALSE, ...) {
  if (is.null(spatial_object)) {
    stop("Parameter 'spatial_object' is missing or null, but is required.")
  }

  # if got the methods name and not the token
  if (method %in% names(deconvolution_methods)) {
    method <- deconvolution_methods[[method]]
  }
  method <- tolower(method)
  check_and_install(method)

  # TODO Type checks for the spatial and single cell object
  # General type checks will be performed here, also matrix + annotation handling
  # Method specific processing steps will be located in the switch statement

  if (!is.null(single_cell_object) && !(cell_type_col %in% names(colData(single_cell_object)))) {
    stop(paste0("Provided col name \"", cell_type_col, "\" can not be found in single_cell_object"))
  }

  # convert data
  if (!is.null(single_cell_object)) {
    single_cell_object <- convert_to_sce(single_cell_object)
  }


  deconv <- switch(method,
    rctd = {
      deconvolute_rctd(single_cell_object, cell_type_col, spatial_object, assay_sc = assay_sc, assay_sp = assay_sp, ...)
    },
    spotlight = {
      deconvolute_spotlight(spe = spatial_object, model = signature, assay_sp = assay_sp)
    },
    card = {
      deconvolute_card(single_cell_object, spatial_object, cell_type_col=cell_type_col, assay_sc = assay_sc, assay_sp = assay_sp, batch_id_col = batch_id_col)
    },
    spatialdwls = {
      deconvolute_spatial_dwls(spe, signature, assay_sp = assay_sp, ...)
    },

    ##############
    # omnideconv #
    ##############
    autogenes = {
      deconvolute_omnideconv(spatial_object, signature, method = "autogenes", single_cell_object, cell_type_col, batch_id_col = batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose)
    },
    bisque = {
      deconvolute_omnideconv(spatial_object, signature, method = "bisque", single_cell_object, cell_type_col, batch_id_col = batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose)
    },
    bseqsc = {
      deconvolute_omnideconv(spatial_object, signature, method = "bseqsc", single_cell_object, cell_type_col, batch_id_col = batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose)
    },
    cibersortx = {
      deconvolute_omnideconv(spatial_object, signature, method = "cibersortx", single_cell_object, cell_type_col, batch_id_col = batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose)
    },
    cdseq = {
      deconvolute_omnideconv(spatial_object, signature, method = "cdseq", single_cell_object, cell_type_col, batch_id_col = batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose)
    },
    cpm = {
      deconvolute_omnideconv(spatial_object, signature, method = "cpm", single_cell_object, cell_type_col, batch_id_col = batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose)
    },
    dwls = {
      deconvolute_omnideconv(spatial_object, signature, method = "dwls", single_cell_object, cell_type_col, batch_id_col = batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose)
    },
    momf = {
      deconvolute_omnideconv(spatial_object, signature, method = "momf", single_cell_object, cell_type_col, batch_id_col = batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose)
    },
    music = {
      deconvolute_omnideconv(spatial_object, signature, method = "music", single_cell_object, cell_type_col, batch_id_col = batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose)
    },
    scaden = {
      deconvolute_omnideconv(spatial_object, signature, method = "scaden", single_cell_object, cell_type_col, batch_id_col = batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose)
    },
    scdc = {
      deconvolute_omnideconv(spatial_object, signature, method = "scdc", single_cell_object, cell_type_col, batch_id_col = batch_id_col, assay_sc = assay_sc, assay_sp = assay_sp, verbose = verbose)
    },

    ################
    # immunedeconv #
    ################
    mcp_counter = {
      deconvolute_immunedeconv(spatial_object, method = "mcp_counter", assay_sp = assay_sp,  ...)
    },
    epic = {
      deconvolute_immunedeconv(spatial_object, method = "epic", assay_sp = assay_sp, ...)
    },
    quantiseq = {
      deconvolute_immunedeconv(spatial_object, method = "quantiseq", assay_sp = assay_sp, ...)
    },
    xcell = {
      deconvolute_immunedeconv(spatial_object, method = "xcell", assay_sp = assay_sp, ...)
    },
    cibersort = {
      deconvolute_immunedeconv(spatial_object, method = "cibersort", assay_sp = assay_sp, ...)
    },
    cibersort_abs = {
      deconvolute_immunedeconv(spatial_object, method = "cibersort_abs", assay_sp = assay_sp, ...)
    },
    timer = {
      deconvolute_immunedeconv(spatial_object, method = "timer", assay_sp = assay_sp, ...)
    },
    consensus_tme = {
      deconvolute_immunedeconv(spatial_object, method = "consensus_tme", assay_sp = assay_sp, ...)
    },
    abis = {
      deconvolute_immunedeconv(spatial_object, method = "abis", assay_sp = assay_sp, ...)
    },
    estimate = {
      deconvolute_immunedeconv(spatial_object, method = "estimate", assay_sp = assay_sp, ...)
    },

    ######################
    # immunedeconv mouse #
    ######################
    mmcp_counter = {
      deconvolute_immunedeconv_mouse(spatial_object, method = "mmcp_counter", assay_sp = assay_sp, ...)
    },
    seqimmucc = {
      deconvolute_immunedeconv_mouse(spatial_object, method = "seqimmucc", assay_sp = assay_sp, ...)
    },
    dcq = {
      deconvolute_immunedeconv_mouse(spatial_object, method = "dcq", assay_sp = assay_sp, ...)
    },
    base = {
      deconvolute_immunedeconv_mouse(spatial_object, method = "base", assay_sp = assay_sp, ...)
    }
  )

  message("finished deconvolution")

  # save to object or return table
  if (return_object) {
    return(addResultToObject(spatial_object, deconv))
  } else {
    return(deconv)
  }
}

#' Build Model and Deconvolute in one step
#'
#' @param single_cell_obj Single Cell Object containing reference data to build the model
#' @param spatial_obj SpatialExperiment to be deconvoluted
#' @param method deconvolution method
#' @param cell_type_col column of single_cell_obj containing cell type information
#' @param batch_id_col column of SpatialObject containing batch_id information
#' @param assay_sc the assay of the single cell object to use, default = "counts"
#' @param assay_sp the assay of the spatialExperiment to use, default = "counts"
#' @param return_object if true return anotation spatialExperiment, if false return table
#' @param verbose output more information
#' @param ... further parameters passed to the selected function
#' @export
build_and_deconvolute <- function(single_cell_obj, spatial_obj, method = NULL, cell_type_col = "cell_ontology_class", batch_id_col = NULL, assay_sc = "counts", assay_sp = "counts", return_object = TRUE, verbose = FALSE, ...) {

  # TODO useful checks

  signature <- build_model(
    single_cell_obj,
    cell_type_col = cell_type_col,
    spatial_object = spatial_obj,
    method = method,
    batch_id_col = batch_id_col,
    assay_sc = assay_sc,
    assay_sp = assay_sp,
    verbose = verbose,
    ...
  )

  deconv <- deconvolute(
    spatial_object = spatial_obj,
    signature = signature,
    single_cell_object = single_cell_obj,
    cell_type_col = cell_type_col,
    batch_id_col = batch_id_col,
    assay_sc = assay_sc,
    assay_sp = assay_sp,
    verbose = verbose,
    return_object = return_object,
    ...
  )
  return(deconv)
}
