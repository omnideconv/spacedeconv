#' Build Model Spatial DWLS
#' @param sce SingleCellExperiment
#' @param assay_sc Single Cell Object assay to use
#' @param marker_method provide list of marker genes or method to calculate markers (scran, gini, mast)
#' @param cell_type_col column of sce containing cell type information
#' @param dim_method dimension reduction method
#' @param cluster_method cluster method to  use when calculating marker genes
#' @param ... additional paramters
build_model_spatial_dwls <- function(sce, assay_sc = "counts", marker_method = "scran", cell_type_col = "cell_ontology_class", dim_method = "pca", cluster_method="leiden", ...) {


  # TODO Checks
  if (!checkCol(sce, cell_type_col)) {
    stop("cell_type_col not available")
  }

  # spatial expression
  scExpression <- SummarizedExperiment::assay(sce, assay_sc) %>% as("dgCMatrix")


  # markers
  # NOTE the makeSignMatrix function requires a list of genes to use in the signature
  # They just mention "e.g. markers", so i thought to use the marker gene function provided by giotto
  if (marker_method %in% c("scran", "gini", "mast")) {
    message("Calculating markers by using the method: ", marker_method)
    # turn into giotto for markers!!
    obj <- Giotto::createGiottoObject(scExpression)

    obj <- doGiottoWorkflow(obj, dim_method = dim_method, cluster_method = cluster_method)


    message("Calculating Markers")
    # calculate based on selection
    # TODO add assay to use or select default one by handling all that during Giotto object creation
    # TODO Cluster, probably just add _clus
    markers = Giotto::findMarkers(obj,  method = marker_method, cluster_column = cluster_method)

  } else {
    message("Using the provided marker genes")
    markers <- marker_method
  }

  # get cell type vector from object
  cell_types <- as.vector(sce[[cell_type_col]])

  # TODO Check if markers are necessary!!!
  signature <- Giotto::makeSignMatrixDWLSfromMatrix(matrix = scExpression, sign_gene = markers, cell_type_vector = cell_types)

  return(signature)
}

#' Deconvolute Spatial DWLS
#'
#' @param spe Spatial Experiment
#' @param signature Signature
#' @param assay_sp Assay of SpatialExperiment to use
#' @param ... additional parameters
deconvolute_spatial_dwls <- function(spe, signature, assay_sp = "counts", ...) {

  # TODO checks

  # create Giotto Object
  spExpression <- SummarizedExperiment::assay(spe, assay_sp) %>% as("dgCMatrix")
  spCoords <- SpatialExperiment::spatialCoords(spe)

  obj <- Giotto::createGiottoObject(raw_exprs = spExpression, spatial_locs = spCoords)

  obj <- doGiottoWorkflow(obj, ...)

  deconvolution <- Giotto::runDWLSDeconv(obj, sign_matrix = signature, n_cell = 10, return_gobject = FALSE)

  return(deconvolution)
}

#' Perform Giotto Workflow of Normalization, HVG, Dimension Reduction and Clustering
#'
#' @param obj Giotto Object
#' @param calculateHVG Wheter to calculate HVG
#' @param dim_method dimension reduction method to use: pca, tsne, umap
#' @param cluster_method cluster method to use c(leiden, kmeans, hclust, louvain)
doGiottoWorkflow <- function(obj, calculateHVG = TRUE, dim_method = "pca", cluster_method = "leiden") {
  if (!class(obj)[[1]] == "giotto") {
    stop("Object needs to be a Giotto object!")
  }

  if (!dim_method %in% c("pca", "tsne", "umap"))

    message("Starting Giotto Workflow")

  # TODO check if all those steps have already been performed???
  obj <- Giotto::normalizeGiotto(obj)

  # HVG
  # NOTE: the workflow fails often when not calculating hvgs, so this should be default!!!!!
  if (calculateHVG) {
    message("Giotto: Calculating HVGs")
    obj <- Giotto::calculateHVG(obj, show_plot = FALSE, return_plot = FALSE, save_plot = FALSE)
  }


  # Dimension Reduction
  message("Giotto: Performing Dimension Reduction")
  obj <- Giotto::runPCA(obj) # run anayway, either pca is used or it needs to be run for tnse and umap?
  obj <- switch(dim_method,
                pca = {obj}, # just return the object, already performed
                tsne = {
                  Giotto::runtSNE(obj) # run based on pca, could also be done directly
                },
                umap = {
                  Giotto::runUMAP(obj) # run based on pca
                }
  )

  # Create Network using the performed dimension reduction
  message("Giotto: Calculating Network")
  obj <- Giotto::createNearestNetwork(
    gobject = obj,
    dim_reduction_to_use = dim_method,
    dim_reduction_name = dim_method
  )


  # TODO add parameters for all methods
  # Clustering
  message("Giotto: Clustering")
  obj <- switch(cluster_method,
                leiden = {
                  Giotto::doLeidenCluster(obj)
                },
                louvain = {
                  Giotto::doLouvainCluster(obj)
                },
                kmeans = {
                  Giotto::doKmeans(obj,
                                   dim_reduction_to_use = dim_method,
                                   dim_reduction_name = dim_method
                  )
                },
                hclust = {
                  Giotto::doHclust(obj,
                                   dim_reduction_to_use = dim_method,
                                   dim_reduction_name = dim_method
                  )
                },
  )

  message("Giotto: Finished")
  return(obj)
}
