#' Print info about dataset
#'
#' @param sce singleCellExperiment
#' @param spe SpatialExperiment
#' @param signature signature matrix
#' @param assay assay to use for computation
#' @export
print_info <- function(sce = NULL, spe = NULL, signature = NULL, assay="counts") {
  # check for correct class

  if (!is.null(sce) & !is(sce, "SingleCellExperiment")) {
    stop("Your single Cell object is not a valid datatype")
  }

  if (!is.null(spe) & !is(spe, "SpatialExperiment")) {
    stop("The spatial object is not a SpatialExperiment")
  }

  # exclamation, rocket, check, cross mark, ok, hourglass, entry
  # sce
  if (!is.null(sce)) {

    if ("counts" %in% assayNames(sce)){
      assaysce <- "counts"
    } else {
      assaysce <- assayNames(sce)[1]
      cli::cli_alert_info(paste("Using assay", assaysce, "for single cell data"))
    }

    cli::cli_h3("Single Cell")
    cli::cli_text("Assays: {.val {assayNames(sce)}}")
    cli::cli_text("Genes: {.val {nrow(sce)}}")

    genes0 <- sum(DelayedArray::rowSums(assay(sce, assaysce)) == 0)
    percentGenes0 <- round(genes0 / nrow(sce) * 100, 2)
    cli::cli_alert("without expression: {.val {genes0}} ({percentGenes0}%)")

    cli::cli_text("Cells: {.val {ncol(sce)}}")
    umi <- colSums(assay(sce, assaysce))
    cells0 <- sum(umi == 0)
    percentCells0 <- round(cells0 / ncol(sce) * 100, 2)
    cli::cli_alert("without expression: {.val {cells0}} ({percentCells0}%)")
    cli::cli_text("Umi count range: {.val {min(umi)}} - {.val {max(umi)}}")


    if (is.null(rownames(sce))) {
      cli::cli_alert_danger("Rownames not set")
    } else {
      cli::cli_alert_success("Rownames set")
    }

    if (is.null(colnames(sce))) {
      cli::cli_alert_danger("Colnames not set")
    } else {
      cli::cli_alert_success("Colnames set")
    }
  }


  # spe

  if (!is.null(spe)) {

    if ("counts" %in% assayNames(spe)){
      assayspe <- "counts"
    } else {
      assayspe <- assayNames(spe)[1]
      cli::cli_alert_info(paste("Using assay", assayspe, "for single cell data"))
    }

    message()
    cli::cli_h3("Spatial")
    cli::cli_text("Assays: {.val {assayNames(spe)}}")
    cli::cli_text("Genes: {.val {nrow(spe)}}")

    genes0 <- sum(DelayedArray::rowSums(assay(spe, assayspe)) == 0)
    percentGenes0 <- round(genes0 / nrow(spe) * 100, 2)
    cli::cli_alert("without expression: {.val {genes0}} ({percentGenes0}%)")

    cli::cli_text("Spots: {.val {ncol(spe)}}")

    if ("in_tissue" %in% names(colData(spe))) {
      spotsInTissue <- length(colData(spe)$in_tissue)
      cli::cli_text("Spots under tissue: {.val {spotsInTissue}} ({round(spotsInTissue/ncol(spe)*100, 2)}%)")
    } else {
      cli::cli_alert_warning("Spatial Object does not contain tissue presence annotation")
    }
    medianGenesPerSpot <- median(DelayedArray::colSums(assay(spe, assayspe) >= 1))
    cli::cli_text("Median Genes Per Spot: {.val {medianGenesPerSpot}}")


    umi <- colSums(assay(spe, assayspe))
    umiBelow500 <- umi[umi < 500]
    umiBelow500Percent <- round(length(umiBelow500) / nrow(spe) * 100, 2)
    cells0 <- sum(umi == 0)
    percentCells0 <- round(cells0 / ncol(spe) * 100, 2)
    cli::cli_alert("without expression: {.val {cells0}} ({percentCells0}%)")
    cli::cli_text("Umi count range: {.val {min(umi)}} - {.val {max(umi)}}")
    cli::cli_text("Spots with UMI  count below 500: {.val {length(umiBelow500)}} ({umiBelow500Percent}%)")


    if (is.null(rownames(spe))) {
      cli::cli_alert_danger("Rownames not set")
    } else {
      cli::cli_alert_success("Rownames set")
    }

    if (is.null(colnames(spe))) {
      cli::cli_alert_danger("Colnames not set")
    } else {
      cli::cli_alert_success("Colnames set")
    }
  }

  # signature
  if (!is.null(signature)) {
    cli::cli_h3("Signature")

    cli::cli_text("Number of Genes: {.val {nrow(signature)}}")
    cli::cli_text("Number of Cell Types: {.val {ncol(signature)}}")
    cli::cli_alert("{.val {colnames(signature)}}")

    if (!is.null(spe)) {
      overlapGenes <- sum(rownames(spe) %in% rownames(signature))
      overlapGenesPercent <- round(overlapGenes / length(rownames(signature)) * 100, 2)
      cli::cli_alert_info("{.val {overlapGenes}} ({overlapGenesPercent}%) signature genes are available in spatial object")
    }
  }
}
