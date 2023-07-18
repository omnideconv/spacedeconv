#' Convert Input Files to SingleCellExperiment
#'
#' @param obj object provided by the user, will be converted to sce
#'
convert_to_sce <- function(obj) {
  sce <- NULL

  if (!is.null(obj)) {
    # check object type and convert if necessary
    if (is(obj, "SingleCellExperiment")) {
      # class(obj)[[1]] == "SingleCellExperiment"
      sce <- obj
    } else if (is(obj, "AnnDataR6")) {
      # class(obj)[[1]] == "AnnDataR6"
      sce <- anndata_to_singlecellexperiment(obj)
      message("Converted Anndata to SCE")
    } else if (is(obj, "Seurat")) {
      # (class(obj)[[1]] == "Seurat"
      sce <- Seurat::as.SingleCellExperiment(obj)
      message("Converted Seurat to SCE")
    }
  }

  return(sce)
}

#' Convert AnnData to SingleCellExperiment
#'
#' @param ad AnnData object
#'
#' @return SingleCellObject
#' Thanks to Omnideconv
#' @export
anndata_to_singlecellexperiment <- function(ad) {
  # anndata_checkload()
  ad <- ad$transpose()
  X_mat <- ad$X
  rownames(X_mat) <- ad$obs_names
  colnames(X_mat) <- ad$var_names


  assays_list <- list()
  assays_list[["X"]] <- X_mat
  assays_list <- c(assays_list, lapply(ad$layers$keys(), function(x) ad$layers[x]))
  names(assays_list) <- c("X", ad$layers$keys())

  meta_list <- list()
  meta_list[ad$uns_keys()] <- ad$uns

  rowdata <- ad$obs
  if (length(ad$obsm) == nrow(rowdata)) {
    rowdata$obsm <- ad$obsm
  }


  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_list,
    rowData = rowdata,
    colData = ad$var,
    reducedDims = ad$varm,
    metadata = meta_list,
    rowPairs = ad$obsp,
    colPairs = ad$varp
  )

  return(sce)
}
#' Convert Seurat to SpatialExperiment
#' @param seurat Seurat Object
#' @returns SpatialExperiment
#'
#' @export
seurat_to_spatialexperiment <- function(seurat) {
  # sce = Seurat::as.SingleCellExperiment(seurat)

  # assay_names = SummarizedExperiment::assays(sce) %>% names()

  images <- Seurat::Images(seurat)

  raster <- SpatialExperiment::SpatialImage(as.raster(seurat@images[[images]]@image))

  img <- S4Vectors::DataFrame(
    sample_id = "sample01", image_id = "lowres", data = I(list(raster)),
    scaleFactor = seurat@images[[images]]@scale.factors[["lowres"]]
  )


  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = Seurat::GetAssayData(seurat, "Spatial", slot = "counts")),
    # spatialCoords = as.matrix(Seurat::GetTissueCoordinates(seurat)) ,
    spatialCoords = as.matrix(seurat@images[[images]]@coordinates[c("imagerow", "imagecol")]),
    imgData = img
  )

  return(spe)
}

# spe_to_seurat <- function(spe){
#   obj <- Seurat::CreateSeuratObject(counts = , assay = "Spatial", project = )
# }

#' Convert AnnData to SpatialExperiment
#' @param ad anndata object
#' @export
anndata_to_spatialexperiment <- function(ad) {
  expr <- Matrix::t(Matrix::Matrix(ad$X, sparse = T))

  sample_metadata <- ad$obs
  gene_metadata <- ad$var

  spatial_coords <- ad$obsm$spatial


  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = expr),
    colData = sample_metadata,
    rowData = gene_metadata,
    spatialCoords = spatial_coords,
    sample_id = as.character(ad$obs$sample)
  )

  images <- names(ad$uns[["spatial"]])
  df <- NULL

  for (sample in levels(ad$obs$sample)) {
    matchingSample <- images[grep(sample, images)]

    lowres <- ad$uns$spatial[[matchingSample]]$images$lowres
    scalefactor <- ad$uns$spatial[[matchingSample]]$scalefactors$tissue_lowres_scalef

    # spot_dim <- ad$uns$spatial[[tmp]]$scalefactors$spot_diameter_fullres

    #plot(as.raster(lowres))
    if (is.null(df)) {
      df <- DataFrame(sample_id = sample, image_id = "lowres", data = I(list(SpatialImage(as.raster(lowres)))), scaleFactor = scalefactor)
    } else {
      df <- rbind(df, data.frame(sample_id = sample, image_id = "lowres", data = I(list(SpatialImage(as.raster(lowres)))), scaleFactor = scalefactor))
    }
  }

  SpatialExperiment::imgData(spe) <- df

  spatialCoordinates <- ad$obsm$spatial

  colnames(spatialCoordinates) <- c("pxl_col_in_fullres", "pxl_row_in_fullres")

  SpatialExperiment::spatialCoords(spe) <- spatialCoordinates


  return(spe)
}


#' Convert Spatial Experiment to AnnData
#' @param spe Spatial experiment
#' @param assay assay to use
#' @export
spe_to_ad <- function(spe, assay = "counts") {
  if (is.null(spe)) {
    stop("Parameter 'spe' is missing or null, but is required")
  }

  if (!assay %in% names(SummarizedExperiment::assays(spe))) {
    stop("Requested assay not available in object")
  }

  ad <- anndata::AnnData(
    X = Matrix::t(SummarizedExperiment::assay(spe, assay)),
    var = as.data.frame(SingleCellExperiment::rowData(spe)),
    obs = as.data.frame(SingleCellExperiment::colData(spe))
  )

  return(ad)
}
