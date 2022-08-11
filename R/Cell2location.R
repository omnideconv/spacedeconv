#' Build Model Cell2location
#'
#' @param sc_obj SingleCellExperiment containing reference expression
#' @param epochs number of epochs to train the network
#' @param cell_count_cutoff cell2location parameter
#' @param cell_percentage_cutoff cell2location parameter
#' @param nonz_mean_cutoff cell2location parameter
build_model_cell2location <- function(sc_obj, epochs = 20, cell_count_cutoff=5, cell_percentage_cutoff=0.03, nonz_mean_cutoff=1.12){


  # build anndata, gene names as rownames


  # source python script
  reticulate::source_python("~/SpaceDeconv/inst/python/cell2location.py")

  model = py_build_model_cell2location(adat_ref,
                                       epochs = epochs,
                                       cell_count_cutoff = cell_count_cutoff,
                                       cell_percentage_cutoff = cell_percentage_cutoff,
                                       nonz_mean_cutoff = nonz_mean_cutoff)

}

#' Deconvolute Cell2location

deconvolute_cell2location <- function(){

}
