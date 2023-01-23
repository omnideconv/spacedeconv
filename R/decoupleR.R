#' Obtain a decoupleR reference
#' @param method method to use, progeny or dorothea
#' @param organism which organism
#' @param n_genes number genes to return, for progeny
#' @export
get_decoupleR_reference <- function(method="progeny", organism ="human", n_genes=500){
  if (method=="progeny"){
    reference <- decoupleR::get_progeny(organism=organism, top=n_genes)
  } else if (method=="dorothea"){
    reference <- decoupleR::get_dorothea(organism = organism) ############## missing parameters!
  } else {
    reference <- NULL
    cli::cli_alert_danger("DecoupleR method not supported")
  }

  return (reference)
}

#' Compute Pathway activities with decoupleR
#' @param spe SpatialExperiment
#' @param reference decoupleR reference
#' @param method calculation method to use
#' @param assay which assay to use
#' @export
compute_decoupleR_activites <- function(spe, reference, method="wsum", assay="counts"){

  cli::cli_progress_step("testing parameter", msg_done = "parameter OK")

  ##### additional checks

  # extract expression matrix
  mat <- assay(spe, assay)


  cli::cli_progress_step(msg="Running decoupleR", msg_done = "Finished")

  if (method=="wsum"){
    res <- run_wsum(mat = mat, network = reference, .source = "source", .target = "target", .mor = "weight", times=100, minsize = 5)
  } else {
    res <- NULL
    cli::cli_alert_danger("DecoupleR method not supported")
  }

  # add results to spe

  df <- res[res$statistic=="corr_wsum", ] #############add parameter

  df$source <- make.names(df$source) ######## okay?

  df <- tidyr::pivot_wider(df, names_from = "source", values_from = "score", id_cols = "condition")

  df$condition <- NULL # remove spot

  df <- attachToken(df, token = "decoupleR")

  colData(spe) <- cbind(colData(spe), df)

  cli::cli_progress_done()

  return (spe)
}
