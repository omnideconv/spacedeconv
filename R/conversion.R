#' Convert Input Files to SingleCellExperiment
#'
#' @param obj object provided by the user, will be converted to sce
convert_to_sce <- function (obj){
  sce = NULL

  if (!is.null(obj)){
    # check object type and convert if necessary
    if (class(obj)[[1]] == "SingleCellExperiment"){
      sce = obj
    } else if(class(obj)[[1]] == "AnnDataR6"){
      sce = anndata_to_singlecellexperiment(obj)
      message("Converted Anndata to SCE")
    } else if (class(obj)[[1]] == "Seurat"){
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
#'
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
