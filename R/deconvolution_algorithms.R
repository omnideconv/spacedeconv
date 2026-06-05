#' Supported Deconvolution Methods
#'
#' Named character vector of supported methods. Names are the display labels and
#' values are the internal method tokens used by the API.
#'
#' @details
#' Second-generation spatial methods: `RCTD`, `SPOTlight`, `CARD`, `spatialDWLS`,
#' `cell2location`, `DOT`, `Rectangle`.\cr
#' First-generation immunedeconv methods: `MCPcounter`, `EPIC`, `quanTIseq`, `xCell`,
#' `CIBERSORT`, `CIBERSORT (abs.)`, `TIMER`, `ConsensusTME`, `ABIS`, `ESTIMATE`.\cr
#' First-generation immunedeconv mouse methods: `mMCPcounter`, `seqImmuCC`, `DCQ`, `BASE`.
#'
#' @export
deconvolution_methods <- c(
  # spatial
  "RCTD" = "rctd",
  "SPOTlight" = "spotlight",
  "CARD" = "card",
  "spatialDWLS" = "spatialdwls",
  "cell2location" = "cell2location",
  "DOT" = "dot",
  "Rectangle" = "rectangle",
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


#' Build a Reference Signature
#'
#' Computes a cell-type signature from annotated scRNA-seq data or dispatches to
#' a method-specific builder. The resulting signature can be passed to
#' `deconvolute()` as the `signature` argument. Some methods (e.g., CARD, DOT,
#' immunedeconv) do not produce a signature and return `NULL`.
#'
#' @param single_cell_obj Single-cell object with cell type annotations.
#' @param cell_type_col Column name containing cell type labels.
#' @param method Signature method; one of `spacedeconv::deconvolution_methods`.
#' @param verbose Print extra progress information.
#' @param spatial_obj SpatialExperiment; required for SPOTlight.
#' @param batch_id_col Batch ID column in the single-cell object (used by some methods).
#' @param assay_sc Single-cell assay to use (default: "counts").
#' @param assay_sp Spatial assay to use (default: "counts").
#' @param ... Additional parameters passed to the selected method.
#'
#' @return A cell-type signature matrix, or `NULL` for methods that build internally.
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
    rectangle = {
      build_model_rectangle()
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
#' Runs the selected deconvolution method on a `SpatialExperiment`. Methods may
#' use a precomputed `signature` (from `build_model()`) or build internally.
#' See `spacedeconv::deconvolution_methods` for the full list of methods.
#'
#' @param spatial_obj A `SpatialExperiment` to deconvolute.
#' @param signature Gene expression signature matrix (if required by the method).
#' @param single_cell_obj A `SingleCellExperiment` with reference cells (if required).
#' @param cell_type_col Column name in `single_cell_obj` containing cell types.
#' @param method Deconvolution method; one of `spacedeconv::deconvolution_methods`.
#' @param batch_id_col Batch ID column in `single_cell_obj` (used by some methods).
#' @param assay_sc Single-cell assay to use (default: "counts").
#' @param assay_sp Spatial assay to use (default: "counts").
#' @param return_object If `TRUE`, return the input `SpatialExperiment` with
#' results added to `colData`; otherwise return the result matrix.
#' @param verbose Print extra progress information.
#' @param ... Further parameters passed to the selected method.
#' @return Deconvolution results as a matrix or a `SpatialExperiment`.
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
    rectangle = {
      deconvolute_rectangle(single_cell_obj = single_cell_obj, spatial_obj = spatial_obj, cell_type_col = cell_type_col, assay_sc = assay_sc, assay_sp = assay_sp, ...)
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

#' Build Model and Deconvolute in One Step
#'
#' Convenience wrapper that builds a signature with `build_model()` and then
#' calls `deconvolute()` using the same inputs.
#'
#' @param single_cell_obj Single-cell object used to build the signature.
#' @param spatial_obj `SpatialExperiment` to deconvolute.
#' @param method Deconvolution method; one of `spacedeconv::deconvolution_methods`.
#' @param cell_type_col Column in `single_cell_obj` with cell type labels.
#' @param batch_id_col Batch ID column (used by some methods).
#' @param assay_sc Single-cell assay to use (default: "counts").
#' @param assay_sp Spatial assay to use (default: "counts").
#' @param return_object If `TRUE`, return an annotated `SpatialExperiment`;
#' otherwise return the result matrix.
#' @param verbose Print extra progress information.
#' @param ... Additional parameters passed to the selected methods.
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
