#' Convert Input Files to SingleCellExperiment
#'
#' @param obj object provided by the user, will be converted to sce
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
  anndata_checkload()
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

#' Convert SpatialExperiment to AnnData
anndata_to_spatialexperiment <- function() {

}



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
