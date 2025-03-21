#' Build Model Cell2location
#'
#' @param single_cell_obj SingleCellExperiment containing reference expression
#' @param epochs number of epochs to train the network
#' @param assay_sc assay to extract
#' @param batch_id_col batch_id_col column
#' @param cell_type_col column of single_cell_obj containing cell type labels
#' @param cell_count_cutoff cell2location parameter
#' @param cell_percentage_cutoff cell2location parameter
#' @param nonz_mean_cutoff cell2location parameter
#' @param gpu whether to train on GPU
build_model_cell2location <- function(single_cell_obj, epochs = 250, assay_sc = "counts", batch_id_col = "sample_id", cell_type_col = "celltype_major", cell_count_cutoff = 5, cell_percentage_cutoff = 0.03, nonz_mean_cutoff = 1.12, gpu = TRUE) {
  # build anndata, gene names as rownames

  # init_python()

  ad <- spe_to_ad(single_cell_obj, assay = assay_sc) # using the spatial function

  # source python script
  reticulate::source_python(system.file("python", "cell2location_spacedeconv.py", package = "spacedeconv"))

  model <- py_build_model_cell2location(ad,
    epochs = as.integer(epochs), # int!
    batch_id_col = batch_id_col,
    cell_type_column = cell_type_col,
    cell_count_cutoff = as.integer(cell_count_cutoff), # int!
    cell_percentage_cutoff = cell_percentage_cutoff,
    nonz_mean_cutoff = nonz_mean_cutoff,
    gpu = gpu
  )
}

#' Deconvolute Cell2location
#'
#' @param spatial_obj SpatialExperiment
#' @param signature signature
#' @param epochs training epochs for model
#' @param n_cell cell2location hyperparameter
#' @param alpha cell2location hyperparameter
#' @param gpu whether to use nvidia gpu for training
#' @param result_name token to identify deconvolution results in object, default = "card"
#' @param values relative or absolute, default: relative
deconvolute_cell2location <- function(spatial_obj, signature = NULL, epochs = 30000, n_cell = 10, alpha = 20, gpu = TRUE, result_name = "cell2location", values = "relative") {
  # init_python()

  # TURN INTO ANNDATA
  ad <- spe_to_ad(spatial_obj)

  # source python script
  reticulate::source_python(system.file("python", "cell2location_spacedeconv.py", package = "spacedeconv")) # ("~/spacedeconv/inst/python/cell2location.py")

  deconv <- py_deconvolute_cell2location(
    sp_obj = ad,
    signature = signature, # must be pandas
    epochs = as.integer(epochs),
    n_cell = as.integer(n_cell), alpha = as.integer(alpha), gpu = gpu
  )

  deconv <- attachToken(deconv, result_name)

  if (values == "relative") {
    # abundance_per_spot = rowSums(data.frame(colData(deconv)[, available_results(deconv, method="c2l")]))
    abundance_per_spot <- rowSums(deconv)

    cli::cli_progress_step("Rescaling Cell2location results to relative fractions", msg_done = "Rescaled Cell2location results to relative fractions")
    # for (result in available_results(deconv, method = "c2l")){
    #   colData(deconv)[, result] <- colData(deconv)[, result]/abundance_per_spot
    # }

    for (i in 1:ncol(deconv)) {
      deconv[, i] <- deconv[, i] / abundance_per_spot
    }

    cli::cli_progress_done()
  }

  return(deconv)
}
