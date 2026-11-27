#' Fetch decoupleR Reference Networks
#'
#' Retrieves pathway or transcription factor target networks from `decoupleR` (Progeny, DoRothEA,
#' or CollecTRI).
#'
#' @param method Reference to retrieve: `"progeny"`, `"dorothea"`, or `"collectri"`.
#' @param organism Target organism (e.g., `"human"` or `"mouse"`).
#' @param n_genes Number of top genes for Progeny (only used with `"progeny"`).
#' @param confidence Confidence levels (A-E) to filter DoRothEA targets (only
#' used with `"dorothea"`).
#' @param ... Additional parameters passed to `decoupleR`.
#'
#' @return A data frame containing the requested reference.
#'
#' @export
get_decoupleR_reference <- function(method = "progeny", organism = "human", n_genes = 500, confidence = NULL, ...) {
  requireNamespace("decoupleR")
  cli::cli_rule(left = "spacedeconv")

  cli::cli_progress_step("Getting decoupleR reference", msg_done = "Got decoupleR reference")

  if (method == "collectri" && !is.null(confidence)) {
    warning("Collectri does not use confidence levels, will not use provided confidence vector")
  }

  if (method == "progeny") {
    reference <- fetch_activity_reference("progeny", organism = organism, top = n_genes, ...)
  } else if (method == "dorothea") {
    reference <- fetch_activity_reference("dorothea", organism = organism, ...) ############## missing parameters!
  } else if (method == "collectri") {
    reference <- fetch_activity_reference("collectri", organism = organism, ...)
  } else {
    reference <- NULL
    cli::cli_alert_danger("DecoupleR method not supported")
  }

  cli::cli_progress_done()

  # only use confidence when requesting dorothea
  if (!is.null(confidence) && method == "dorothea") {
    reference <- reference[reference$confidence %in% confidence, ]
  }

  return(reference)
}

#' Compute Pathway and Transcription Factor Activities with decoupleR
#'
#' Computes pathway or transcription factor activities from a `SpatialExperiment`
#' using methods provided by the `decoupleR` package.
#'
#' @param spe `SpatialExperiment` containing expression data.
#' @param reference A data frame with `source` (gene set or regulator), `target`
#' (gene identifier), and `weight` or `mor` (weight or regulatory direction).
#' Custom tables can be supplied directly; built-in references are retrieved with
#' `get_decoupleR_reference()`. If both weight columns are present, `mor` is used.
#' The activity scoring methods AUCell, FGSEA, GSVA, and ORA do not use weights. Gene identifiers must match
#' the expression assay row names and organism. At least five matching genes
#' per set are required. Each source-target pair should occur only once.
#' Custom tables with `weight` and no database-specific metadata use the
#' `decoupleR_` column prefix.
#' @param method Activity method (e.g., `"wmean"`, `"viper"`, `"gsva"`); see
#' `decoupleR` docs for options.
#' @param assay Assay to use for calculations.
#' @param statistic If multiple statistics are produced, select one to store.
#' @param ... Further arguments passed to the selected method.
#'
#' @return `SpatialExperiment` with activity scores added to `colData`.
#' @examples
#' # Illustrative genes only; use a biologically justified reference for analysis.
#' ref <- data.frame(
#'   source = rep("my_signature", 5),
#'   target = c("STAT1", "IRF1", "CXCL9", "CXCL10", "IDO1"),
#'   weight = rep(1, 5)
#' )
#' # Use known regulatory weights instead of 1 when available.
#' \dontrun{
#' # spe must contain a normalized expression assay named "cpm".
#' spe <- compute_activities(spe, reference = ref, method = "ulm", assay = "cpm")
#' }
#' @export
compute_activities <- function(spe, reference, method = "wmean", assay = "cpm", statistic = NULL, ...) {
  requireNamespace("decoupleR")
  cli::cli_rule(left = "spacedeconv")


  cli::cli_progress_step("testing parameter", msg_done = "parameter OK")

  ##### additional checks

  # extract expression matrix
  mat <- assay(spe, assay)

  # dorothea: mor
  # progeny: weight

  if ("mor" %in% colnames(reference)) {
    values <- "mor"
  } else if ("weight" %in% colnames(reference)) {
    values <- "weight"
  } else {
    cli::cli_alert_danger("Reference not applicable")
    stop()
  }

  cli::cli_progress_step(msg = "Running decoupleR", msg_done = "Finished")

  if (method == "wsum") {
    res <- run_wsum(mat = mat, network = reference, .source = "source", .target = "target", .mor = values, times = 100, minsize = 5)
  } else if (method == "aucell") {
    res <- run_aucell(mat = as.matrix(mat), network = reference, .source = "source", .target = "target", nproc = 4, seed = 42, minsize = 5) # missing params
  } else if (method == "fgsea") {
    res <- run_fgsea(mat = mat, network = reference, .source = "source", .target = "target", times = 100, nproc = 4, seed = 42, minsize = 5)
  } else if (method == "gsva") {
    res <- run_gsva(mat = mat, network = reference, .source = "source", .target = "target", minsize = 5) # missing parameters
  } else if (method == "mdt") {
    res <- run_mdt(mat = mat, network = reference, .source = "source", .target = "target", .mor = values, minsize = 5) # missing parameters
  } else if (method == "mlm") {
    res <- run_mlm(mat = mat, network = reference, .source = "source", .target = "target", .mor = values, minsize = 5) # missing parameters
  } else if (method == "ora") {
    res <- run_ora(mat = mat, network = reference, .source = "source", .target = "target", minsize = 5) # missing parameters
  } else if (method == "udt") {
    res <- run_udt(mat = mat, network = reference, .source = "source", .target = "target", .mor = values, minsize = 5) # missing parameters
  } else if (method == "ulm") {
    res <- run_ulm(mat = mat, network = reference, .source = "source", .target = "target", .mor = values, minsize = 5) # missing parameters
  } else if (method == "viper") {
    res <- run_viper(mat = mat, network = reference, .source = "source", .target = "target", .mor = values, minsize = 5) # missing parameters
  } else if (method == "wmean") {
    res <- run_wmean(mat = mat, network = reference, .source = "source", .target = "target", .mor = values, minsize = 5) # missing parameters
  } else {
    res <- NULL
    cli::cli_alert_danger("DecoupleR method not supported")
  }

  # add results to spe

  # df <- res[res$statistic=="corr_wsum", ] #############add parameter

  tmp <- unique(res$statistic)
  if (length(tmp) > 1) {
    if (!is.null(statistic)) {
      tmp <- statistic
    }
    cli::cli_alert(paste("Calculated multiple results, using", tmp[[1]], " , other available results: ", tmp[-1]))
    df <- res[res$statistic == tmp[[1]], ]
  } else {
    df <- res
  }

  df$source <- make.names(df$source) ######## okay?

  df <- tidyr::pivot_wider(df, names_from = "source", values_from = "score", id_cols = "condition")

  df <- df[match(colnames(spe), df$condition), , drop = FALSE]
  df$condition <- NULL # IDs have been aligned to the input spots

  # check if dorothea or progeny reference, for naming
  if ("p_value" %in% names(reference)) {
    decouple_tool <- "progeny"
  } else if ("confidence" %in% names(reference)) {
    decouple_tool <- "dorothea"
  } else if ("mor" %in% names(reference)) {
    decouple_tool <- "collectri"
  } else {
    decouple_tool <- "decoupleR"
  }

  df <- attachToken(df, token = decouple_tool)

  colData(spe) <- cbind(colData(spe), df)

  cli::cli_progress_done()

  return(spe)
}

# Isolate the remote-data boundary so tests do not depend on a live service.
fetch_activity_reference <- function(method, organism, ...) {
  switch(method,
    progeny = decoupleR::get_progeny(organism = organism, ...),
    dorothea = decoupleR::get_dorothea(organism = organism, ...),
    collectri = decoupleR::get_collectri(organism = organism, ...))
}
