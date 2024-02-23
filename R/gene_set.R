#' Calculate Gene Set Score
#'
#' @param spe SpatialExperiment object
#' @param genes vectors of genes to calculate gene set score
#' @param assay name of the assay to use, default is "cpm"
#' @param name name of the result column
#'
#' @returns updated SpatialExperiment object including the gene set score
#' @export
gene_set_score <- function(spe, genes = NULL, assay = "cpm", name = "geneSet") {
  # Ensure spe is a SpatialExperiment object
  if (!inherits(spe, "SpatialExperiment")) {
    stop("spe must be a SpatialExperiment object.")
  }

  # Check if genes is provided and is a vector of gene names
  if (is.null(genes) || !is.vector(genes) || !all(is.character(genes))) {
    stop("genes must be a non-null vector of gene names.")
  }

  # Ensure the name is a non-empty string
  if (!is.character(name) || length(name) != 1) {
    stop("name must be a non-empty single string.")
  }

  # Subset the SpatialExperiment object by genes
  spe_subset <- spe[rownames(spe) %in% genes, ]

  # Check if the subset is empty
  if (nrow(spe_subset) == 0) {
    stop("No gene available in Spatial Data.")
  }

  # Check if the specified assay is available, fallback to counts if not
  if (assay %in% assayNames(spe)) {
    tmp <- assay(spe_subset, assay)
  } else {
    cli::cli_alert_warning("Requested Assay not in object, using counts instead")
    tmp <- assay(spe_subset, "counts")
  }

  # Calculate the gene set score
  results <- colSums(log(tmp + 1))

  # Add the results as a new column to colData
  colData(spe)[[name]] <- results

  # Return the updated SpatialExperiment object
  return(spe)
}
