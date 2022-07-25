  #' Build Model Spatial DWLS
  #' @param sce SingleCellExperiment
  #' @param assay_sc Single Cell Object assay to use
  #' @param markers provide list of marker genes or method to calculate markers (default, scran, gini, mast)
  #' @param cell_type_col column of sce containing cell type information
  #' @param ... additional paramters
  build_model_spatial_dwls <- function(sce, assay_sc = "counts", marker_method="default", cell_type_col = "cell_ontology_class", ...){


    # TODO Checks
    if (!(cell_type_col) %in% names(SingleCellExperiment::colData(sce))){
      stop("cell_type_col not available")
    }

    # spatial expression
    scExpression <- SummarizedExperiment::assay(sce, assay_sc) %>% as("dgCMatrix")



    # markers
    if (marker_method %in% c("default", "scran", "gini", "mast")){
      message("Calculating markers by using the method: ", marker_method)
      # turn into giotto for markers!!
      obj <- Giotto::createGiottoObject(scExpression)

      obj <- doGiottoWorkflow(obj)


      message("Calculating Markers")
      # calculate based on selection
      if (marker_method=="default"){

        markers = Giotto::findMarkers(obj, cluster_column = "leiden_clus")
      }
    } else {
      message("Using the provided marker genes")
      markers = marker_method
    }



    cell_types <- as.vector(sce[[cell_type_col]]) # get cell types

    signature <- Giotto::makeSignMatrixDWLSfromMatrix(matrix = scExpression, sign_gene = markers, cell_type_vector = cell_types)
  }

  #' Deconvolute Spatial DWLS
  #'
  #' @param spe Spatial Experiment
  #' @param signature Signature
  #' @param assay_sp Assay of SpatialExperiment to use
  #' @param ... additional parameters
  deconvolute_spatial_dwls <- function(spe, signature, assay_sp = "counts", ...){

    # create Giotto Object
    spExpression <- SummarizedExperiment::assay(spe, assay_sp) %>% as("dgCMatrix")
    spCoords <- SpatialExperiment::spatialCoords(spe)

    obj <- Giotto::createGiottoObject(raw_exprs = spExpression, spatial_locs = spCoords)

    obj <- doGiottoWorkflow(obj, ...)

    deconvolution <- Giotto::runDWLSDeconv(obj, sign_matrix = signature, n_cell = 10, return_gobject = FALSE)

    return (deconvolution)
  }

  #' Perform Giotto Workflow of Normalization, HVG, Dimension Reduction and Clustering
  #'
  #' @param obj Giotto Object
  #' @param calculateHVG Wheter to calculate HVG
  #' @param dim_method dimension reduction method to use, pca, tsne, umap
  #' @param cluster_method cluster method to use c(leiden, kmeans, hclust, louvain)

  doGiottoWorkflow <- function(obj, calculateHVG = TRUE, dim_method = "pca", cluster_method = "leiden"){
    if (!class(obj)[[1]] =="giotto"){
      stop("Object needs to be a Giotto object!")
    }

    obj <- Giotto::normalizeGiotto(obj)

    # HVG
    if (calculateHVG){
      message("Calculating HVGs")
      obj <- Giotto::calculateHVG(obj, show_plot = FALSE, return_plot = FALSE, save_plot = FALSE)
    }


    # Dimension Reduction
    message("Performing Dimension Reduction")
    obj <- switch(dim_method,
                  pca = {Giotto::runPCA(obj)},
                  tsne= {},
                  umap = {})

    # Create Network
    message("Calculating Network")
    obj <- Giotto::createNearestNetwork(obj)

    # Clustering
    message("Clustering")
    obj <- Giotto::doLeidenCluster(obj)

    message("Finished")
    return (obj)

  }
