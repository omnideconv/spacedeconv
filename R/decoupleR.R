#' Obtain a decoupleR reference
#' @param method method to use, progeny, dorothea or collectri
#' @param organism which organism
#' @param n_genes number genes to return, for progeny
#' @param confidence condfidence level for transcription factor reference, vector of levels to include
#' @param ... additional parameters to pass to the methods
#' @export
get_decoupleR_reference <- function(method = "progeny", organism = "human", n_genes = 500, confidence = NULL) {
  requireNamespace("decoupleR")
  cli::cli_rule(left = "spacedeconv")

  cli::cli_progress_step("Getting decoupleR reference", msg_done = "Got decoupleR reference")

  if (method == "progeny") {
    reference <- decoupleR::get_progeny(organism = organism, top = n_genes, ...)
  } else if (method == "dorothea") {
    reference <- decoupleR::get_dorothea(organism = organism, ...) ############## missing parameters!
  } else if (method=="collectri"){
    reference <- decoupleR::get_collectri(organism = organism, ...)
  } else {
    reference <- NULL
    cli::cli_alert_danger("DecoupleR method not supported")
  }

  cli::cli_progress_done()

  if (!is.null(confidence)) {
    reference <- reference[reference$confidence %in% confidence, ]
  }

  return(reference)
}

#' Compute Pathway activities with decoupleR
#' @param spe SpatialExperiment
#' @param reference decoupleR reference
#' @param method calculation method to use
#' @param assay which assay to use
#' @param statistic select a sub results in case methods produce mutliple ones
#' @param ... further arguments passed to the methods
#' @export
compute_decoupleR_activities <- function(spe, reference, method = "wmean", assay = "cpm", statistic = NULL, ...) {
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

  df$condition <- NULL # remove spot

  # check if dorothea or progeny reference, for naming
  if ("p_value" %in% names(reference)) {
    decouple_tool <- "progeny"
  } else if ("mor" %in% names(reference)) {
    decouple_tool <- "dorothea"
  } else {
    decouple_tool <- "decoupleR"
  }

  df <- attachToken(df, token = decouple_tool)

  colData(spe) <- cbind(colData(spe), df)

  cli::cli_progress_done()

  return(spe)
}
