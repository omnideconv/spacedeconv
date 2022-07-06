#' List of supported deconvolution methods
#' @export
#'
deconvolution_methods <- c(
  "RCTD" = "rctd",
  "SPOTlight" = "spotlight",
  "CARD" = "card",
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
#' @export
build_model <- function(single_cell_object, cell_type_col = "cell_ontology_class", method = NULL, verbose = FALSE) {
  if (is.null(single_cell_object)) {
    stop("Patameter 'single_cell_object' missing or null, but is required")
  }

  # if got the methods name and not the token
  if (method %in% names(deconvolution_methods)) {
    method <- deconvolution_methods[[method]]
  }
  method <- tolower(method)

  # check if cell_type_col in names(colData(sc)^)


  signature <- switch(method,
    rctd = {build_model_rctd()},

    ##############
    # omnideconv #
    ##############
    autogenes = {
      omnideconv::build_model(single_cell_object, cell_type_column_name = cell_type_col, method = "autogenes", verbose = verbose)
    },
    bisque = {
      omnideconv::build_model(single_cell_object, cell_type_column_name = cell_type_col, method = "bisque", verbose = verbose)
    },
    bseqsc = {
      omnideconv::build_model(single_cell_object, cell_type_column_name = cell_type_col, method = "bseqsc", verbose = verbose)
    },
    cibersortx = {
      omnideconv::build_model(single_cell_object, cell_type_column_name = cell_type_col, method = "dwls", verbose = verbose)
    },
    cdseq = {
      omnideconv::build_model(single_cell_object, cell_type_column_name = cell_type_col, method = "cdseq", verbose = verbose)
    },
    cpm = {
      omnideconv::build_model(single_cell_object, cell_type_column_name = cell_type_col, method = "cpm", verbose = verbose)
    },
    dwls = {
      omnideconv::build_model(single_cell_object, cell_type_column_name = cell_type_col, method = "dwls", verbose = verbose)
    },
    momf = {
      omnideconv::build_model(single_cell_object, cell_type_column_name = cell_type_col, method = "momf", verbose = verbose)
    },
    music = {
      omnideconv::build_model(single_cell_object, cell_type_column_name = cell_type_col, method = "music", verbose = verbose)
    },
    scaden = {
      omnideconv::build_model(single_cell_object, cell_type_column_name = cell_type_col, method = "scaden", verbose = verbose)
    },
    scdc = {
      omnideconv::build_model(single_cell_object, cell_type_column_name = cell_type_col, method = "scdc", verbose = verbose)
    },

    ################
    # immunedeconv #
    ################
    mcp_counter = {},
    epic = {},
    quantiseq = {},
    xCell = {},
    cibersort = {},
    cibersort_abs = {},
    timer = {},
    consensus_tme = {},
    abis = {},

    ######################
    # immunedeconv mouse #
    ######################
    mmcp_counter = {},
    seqimmucc = {},
    dcq = {},
    base = {}
  )

  return(signature)
}


#' Deconvolution
#' @param spatial_object A SpatialExperiment
#' @param signature Gene Expression Signature
#' @param single_cell_object A SingleCellExperiment
#' @param cell_type_col Column name of the single_cell_object where the cell type can be found
#' @param method Deconvolution Method to use, see deconvolution_methods() for a full list of available methods
#' @param return_object Return an Object or result Table, TRUE = Object
#' @param verbose display more detailed information
#' @param ... Further parameters passed to the selected deconvolution method
#' @returns The deconvolution result as a table
#' @export
deconvolute <- function(spatial_object, signature = NULL, single_cell_object, cell_type_col = "cell_ontology_class", method = NULL, return_object = TRUE, verbose = FALSE, ...) {
  if (is.null(spatial_object)) {
    stop("Parameter 'spatial_object' is missing or null, but is required.")
  }
  if (is.null(single_cell_object)) {
    stop("Parameter 'single_cell_object' is missing or null, but is required.")
  }

  # if got the methods name and not the token
  if (method %in% names(deconvolution_methods)) {
    method <- deconvolution_methods[[method]]
  }
  method <- tolower(method)

  # TODO Type checks for the spatial and single cell object
  # General type checks will be performed here, also matrix + annotation handling
  # Method specific processing steps will be located in the switch statement

  if (!(cell_type_col %in% names(colData(spatial_object)))) {
    stop(paste0("Provided col name \"", cell_type_col, "\" can not be found in spatial_object"))
  }

  deconv <- switch(method,
    rctd = {
      deconvolute_rctd(single_cell_object, cell_type_col, spatial_object, ...)
    },

    ##############
    # omnideconv #
    ##############
    autogenes = {
    },
    bisque = {
    },
    bseqsc = {
    },
    cibersortx = {
    },
    cdseq = {
    },
    cpm = {
    },
    dwls = {
      bulk <- counts(spatial_object)
      omnideconv::deconvolute(bulk, signature, method = "dwls", verbose = verbose)
    },
    momf = {
    },
    music = {
    },
    scaden = {
    },
    scdc = {
    },

    ################
    # immunedeconv #
    ################
    mcp_counter = {},
    epic = {},
    quantiseq = {},
    xCell = {},
    cibersort = {},
    cibersort_abs = {},
    timer = {},
    consensus_tme = {},
    abis = {},

    ######################
    # immunedeconv mouse #
    ######################
    mmcp_counter = {},
    seqimmucc = {},
    dcq = {},
    base = {}
  )
  return(deconv)
}
