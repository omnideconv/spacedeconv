#' Build a cell2location Model
#'
#' Wraps the cell2location Python workflow via reticulate to train a reference
#' model from single-cell data. Returns the trained model object.
#'
#' @param single_cell_obj `SingleCellExperiment` containing reference expression.
#' @param epochs Number of training epochs.
#' @param assay_sc Single-cell assay to use.
#' @param batch_id_col Batch ID column name.
#' @param cell_type_col Column with cell type labels.
#' @param cell_count_cutoff cell2location parameter.
#' @param cell_percentage_cutoff cell2location parameter.
#' @param nonz_mean_cutoff cell2location parameter.
#' @param gpu Logical; whether to train on GPU.
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

#' Deconvolute with cell2location
#'
#' Runs the cell2location model (via Python/reticulate) and returns a
#' deconvolution matrix, optionally rescaled to relative values.
#'
#' @param spatial_obj `SpatialExperiment`.
#' @param signature Signature/model input expected by cell2location.
#' @param epochs Training epochs for the model.
#' @param n_cell cell2location hyperparameter.
#' @param alpha cell2location hyperparameter.
#' @param gpu Logical; whether to use GPU for training.
#' @param result_name Prefix used to label result columns (default: "cell2location").
#' @param values `"relative"` to rescale to fractions, or `"absolute"` to keep
#' raw outputs.
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
