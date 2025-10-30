#' List of supported deconvolution methods
#'
#' @details Spatial Algorithms \cr
#' `RCTD`, `SPOTlight`, `CARD`, `spatialDWLS`, `Cell2Location`, `DOT`
#'
#' @details First-generation Methods (Immunedeconv) \cr
#' `MCPcounter`, `EPIC`, `quanTIseq`, `xCell`, `CIBERSORT`, `CIBERSORT (abs.)`,
#' `TIMER`, `ConsensusTME`, `ABIS`, `ESTIMATE`
#'
#' @details First-generation Mouse Methods (Immunedeconv) \cr
#' `mMCPcounter`, `seqImmuCC`, `DCQ`, `BASE`
#'
#' This object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
#'
deconvolution_methods <- c(
  # spatial
  "RCTD" = "rctd",
  "SPOTlight" = "spotlight",
  "CARD" = "card",
  "spatialDWLS" = "spatialdwls",
  "cell2location" = "cell2location",
  "DOT" = "dot",
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
  "ConsensusTME" = "consensus_tme",
  "ABIS" = "abis",
  "ESTIMATE" = "estimate",
  # immunedeconv mouse
  "mMCPcounter" = "mmcp_counter",
  "seqImmuCC" = "seqimmucc",
  "DCQ" = "dcq",
  "BASE" = "base"
)

first_gen <- c(
  # immunedeconv
  "MCPcounter" = "mcp_counter",
  "EPIC" = "epic",
  "quanTIseq" = "quantiseq",
  "xCell" = "xcell",
  "CIBERSORT" = "cibersort",
  "CIBERSORT (abs.)" = "cibersort_abs",
  "TIMER" = "timer",
  "ConsensusTME" = "consensus_tme",
  "ABIS" = "abis",
  "ESTIMATE" = "estimate",
  # immunedeconv mouse
  "mMCPcounter" = "mmcp_counter",
  "seqImmuCC" = "seqimmucc",
  "DCQ" = "dcq",
  "BASE" = "base"
)


#' Build a reference signature
#'
#' @description
#' Build_model calculates a reference signature from annotated scRNA-seq expression data
#'
#' @param single_cell_obj Single-cell Object
#' @param cell_type_col Name of the annotation column containing cell type information
#' @param method Signature calculation Algorithm
#' @param verbose Display more information on console
#' @param spatial_obj SpatialExperiment, required for SPOTlight
#' @param batch_id_col column of singleCellExperiment containing batch ids
#' @param assay_sc assay of single cell object to use
#' @param assay_sp assay of spatial object to use
#' @param ... additional parameters passed to the functions
#'
#' @returns cell-type specific expression signature
#'
#' @export
build_model <- function(single_cell_obj, cell_type_col = "cell_ontology_class", method = NULL, verbose = FALSE, spatial_obj = NULL, batch_id_col = NULL, assay_sc = "counts", assay_sp = "counts", ...) {
  cli::cli_rule(left = "spacedeconv")

  cli::cli_progress_step("testing parameter", msg_done = "parameter OK")

  if (is.null(single_cell_obj)) {
    stop("Parameter 'single_cell_obj' missing or null, but is required")
  }

  # check if method null or not supported
  if (is.null(method) || !(method %in% deconvolution_methods)) {
    stop("Parameter 'method' is null or not supported")
  }

  # if got the methods name and not the token
  if (method %in% names(deconvolution_methods)) {
    method <- deconvolution_methods[[method]]
  }

  # convert to sparse matrices
  if (!is.null(spatial_obj)) {
    spatial_obj <- check_datatype(spatial_obj)
  }

  single_cell_obj <- check_datatype(single_cell_obj)

  # check if method available
  method <- tolower(method)
  check_and_install(method)

  # convert data
  single_cell_obj <- convert_to_sce(single_cell_obj)

  # check if rownames and colnames are set
  if (checkRowColumn(single_cell_obj) || (!is.null(spatial_obj) && checkRowColumn(spatial_obj))) {
    stop("Rownames or colnames not set for single_cell_obj or spatial_obj but need to be available!")
  }

  cli::cli_progress_done()

  print_info(sce = single_cell_obj, spe = spatial_obj)

  cli::cli_progress_step("building signature", msg_done = "finished")

  signature <- switch(method,
    rctd = {
      build_model_rctd()
    },
    spotlight = {
      build_model_spotlight(single_cell_obj, cell_type_col, spatial_obj, assay_sc = assay_sc, assay_sp = assay_sp, ...)
    },
    card = {
      build_model_card()
    },
    spatialdwls = {
      build_model_spatial_dwls(single_cell_obj, assay_sc = assay_sc, marker_method = "scran", cell_type_col = cell_type_col, ...)
    },
    cell2location = {
      build_model_cell2location(single_cell_obj, assay_sc = assay_sc, cell_type_col = cell_type_col, batch_id_col = batch_id_col, ...)
    },
    dot = {
      build_model_dot()
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

  cli::cli_progress_done()

  return(signature)
}


#' Deconvolution with spacedeconv
#'
#' Perform cell type deconvolution with spacedeconv.
#'
#' @param spatial_obj A SpatialExperiment
#' @param signature Gene Expression Signature
#' @param single_cell_obj A SingleCellExperiment
#' @param cell_type_col Column name of the single_cell_obj where the cell type can be found
#' @param method Deconvolution Method to use, see deconvolution_methods() for a full list of available methods
#' @param batch_id_col column name of batch ids in single cell object
#' @param assay_sc which single cell assay to use for computation
#' @param assay_sp which spatial assay to use for computation
#' @param return_object Return an Object or result Table, TRUE = Object
#' @param verbose display more detailed information
#' @param ... Further parameters passed to the selected deconvolution method
#' @returns The deconvolution result as a table
#' @export
deconvolute <- function(spatial_obj, signature = NULL, single_cell_obj = NULL,
                        cell_type_col = "cell_ontology_class", method = NULL,
                        batch_id_col = NULL, assay_sc = "counts",
                        assay_sp = "counts", return_object = TRUE,
                        verbose = FALSE, ...) {
  cli::cli_rule(left = "spacedeconv")

  cli::cli_progress_step("testing parameter", msg_done = "parameter OK")




  if (is.null(spatial_obj)) {
    stop("Parameter 'spatial_obj' is missing or null, but is required.")
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

  if (!is.null(single_cell_obj) && !(cell_type_col %in% names(colData(single_cell_obj)))) {
    stop(paste0("Provided col name \"", cell_type_col, "\" can not be found in single_cell_obj"))
  }

  # convert data
  if (!is.null(single_cell_obj)) {
    single_cell_obj <- convert_to_sce(single_cell_obj)
  }

  # check if rownames and colnames are set
  if (checkRowColumn(spatial_obj)) {
    stop("Rownames or colnames not set for single_cell_obj or spatial_obj but need to be available!")
  }

  if (!is.null(single_cell_obj) && checkRowColumn(single_cell_obj)) {
    stop("Rownames or colnames not set for single_cell_obj or spatial_obj but need to be available!")
  }

  if (method %in% first_gen && checkENSEMBL(rownames(spatial_obj))) {
    stop("You requested a first-generation method and your spatial object contains ENSEBML IDs, please provide HGNC Symbols!")
  }

  # convert to sparse matrices
  spatial_obj <- check_datatype(spatial_obj)

  if (!is.null(single_cell_obj)) {
    single_cell_obj <- check_datatype(single_cell_obj)
  }

  cli::cli_progress_done()

  print_info(sce = single_cell_obj, spe = spatial_obj, signature = signature)

  cli::cli_progress_step("deconvoluting", msg_done = "finished")

  deconv <- switch(method,
    rctd = {
      deconvolute_rctd(single_cell_obj, cell_type_col, spatial_obj, assay_sc = assay_sc, assay_sp = assay_sp, ...)
    },
    spotlight = {
      deconvolute_spotlight(spatial_obj = spatial_obj, model = signature, assay_sp = assay_sp, ...)
    },
    card = {
      deconvolute_card(single_cell_obj, spatial_obj, cell_type_col = cell_type_col, assay_sc = assay_sc, assay_sp = assay_sp, batch_id_col = batch_id_col, ...)
    },
    spatialdwls = {
      deconvolute_spatial_dwls(spatial_obj, signature, assay_sp = assay_sp, ...)
    },
    cell2location = {
      deconvolute_cell2location(spatial_obj, signature, ...)
    },
    dot = {
      deconvolute_dot(single_cell_obj = single_cell_obj, spatial_obj = spatial_obj, cell_type_col = cell_type_col, assay_sc = assay_sc, assay_sp = assay_sp, ...)
    },

    ################
    # immunedeconv #
    ################
    mcp_counter = {
      deconvolute_immunedeconv(spatial_obj, method = "mcp_counter", assay_sp = assay_sp, ...)
    },
    epic = {
      require("EPIC") # quick fix, did not work with requireNamespace()
      deconvolute_immunedeconv(spatial_obj, method = "epic", assay_sp = assay_sp, ...)
    },
    quantiseq = {
      deconvolute_immunedeconv(spatial_obj, method = "quantiseq", assay_sp = assay_sp, ...)
    },
    xcell = {
      require("xCell") # quick fix
      deconvolute_immunedeconv(spatial_obj, method = "xcell", assay_sp = assay_sp, ...)
    },
    cibersort = {
      deconvolute_immunedeconv(spatial_obj, method = "cibersort", assay_sp = assay_sp, ...)
    },
    cibersort_abs = {
      deconvolute_immunedeconv(spatial_obj, method = "cibersort_abs", assay_sp = assay_sp, ...)
    },
    timer = {
      if (!hasArg(indications)) {
        cli::cli_alert_warning("Timer requires a cancer type indications vector")
      }
      deconvolute_immunedeconv(spatial_obj, method = "timer", assay_sp = assay_sp, ...)
    },
    consensus_tme = {
      if (!hasArg(indications)) {
        cli::cli_alert_warning("ConsensusTME requires a cancer type indications vector")
      }
      deconvolute_immunedeconv(spatial_obj, method = "consensus_tme", assay_sp = assay_sp, ...)
    },
    abis = {
      deconvolute_immunedeconv(spatial_obj, method = "abis", assay_sp = assay_sp, ...)
    },
    estimate = {
      deconvolute_immunedeconv(spatial_obj, method = "estimate", assay_sp = assay_sp, ...)
    },

    ######################
    # immunedeconv mouse #
    ######################
    mmcp_counter = {
      deconvolute_immunedeconv_mouse(spatial_obj, method = "mmcp_counter", assay_sp = assay_sp, ...)
    },
    seqimmucc = {
      deconvolute_immunedeconv_mouse(spatial_obj, method = "seqimmucc", assay_sp = assay_sp, ...)
    },
    dcq = {
      deconvolute_immunedeconv_mouse(spatial_obj, method = "dcq", assay_sp = assay_sp, ...)
    },
    base = {
      deconvolute_immunedeconv_mouse(spatial_obj, method = "base", assay_sp = assay_sp, ...)
    }
  )


  # save to object or return table
  if (return_object) {
    return(addResultToObject(spatial_obj, deconv))
  } else {
    return(deconv)
  }
  cli::cli_progress_done()
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
build_and_deconvolute <- function(single_cell_obj, spatial_obj, method = NULL,
                                  cell_type_col = "cell_ontology_class",
                                  batch_id_col = NULL, assay_sc = "counts",
                                  assay_sp = "counts", return_object = TRUE,
                                  verbose = FALSE, ...) {
  # TODO useful checks

  cli::cli_rule(left = "spacedeconv")

  cli::cli_progress_step("testing parameter", msg_done = "parameter OK")

  # check if rownames and colnames are set
  if (checkRowColumn(single_cell_obj) || checkRowColumn(spatial_obj)) {
    stop("Rownames or colnames not set for single_cell_obj or spatial_obj but need to be available!")
  }


  cli::cli_progress_step("building signature", msg_done = "finished")

  signature <- build_model(
    single_cell_obj,
    cell_type_col = cell_type_col,
    spatial_obj = spatial_obj,
    method = method,
    batch_id_col = batch_id_col,
    assay_sc = assay_sc,
    assay_sp = assay_sp,
    verbose = verbose,
    ...
  )


  cli::cli_progress_step("deconvoluting", msg_done = "finished")

  deconv <- deconvolute(
    spatial_obj = spatial_obj,
    signature = signature,
    method = method,
    single_cell_obj = single_cell_obj,
    cell_type_col = cell_type_col,
    batch_id_col = batch_id_col,
    assay_sc = assay_sc,
    assay_sp = assay_sp,
    verbose = verbose,
    return_object = return_object,
    ...
  )

  cli::cli_progress_done()

  return(deconv)
}
