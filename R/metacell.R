#' Clean Genes and Cells with Metacell2
#'
#' @param anndata anndata object containing Single Cell Expression Data
#' @param properly_sampled_min_cell_total exclude all cells with UMI count below this threshold
#' @param properly_sampled_max_cell_total exclude all cells with UMI count above this threshold
#' @param properly_sampled_max_excluded_genes_fraction exclude all cells whose sum of excluded data divided by the toal is more then the threshold
#' @param exclude_genes gene names to exlude
#' @param exclude_gene_patterns gene patterns to exclude
#' @param seed seed for reproducibilty
#'
#' @export
clean_genes_and_cells <- function(anndata, properly_sampled_min_cell_total = 800,
                                  properly_sampled_max_cell_total = 8000,
                                  properly_sampled_max_excluded_genes_fraction = 0.1,
                                  exclude_genes = "", exclude_gene_patterns = "", seed = 123456) {
  reticulate::source_python(system.file("python", "metacells.py", package = "spacedeconv"))

  res <- clean_genes_and_cells(anndata,
    properly_sampled_min_cell_total = as.integer(properly_sampled_min_cell_total),
    properly_sampled_max_cell_total = as.integer(properly_sampled_max_cell_total),
    properly_sampled_max_excluded_genes_fraction = properly_sampled_max_excluded_genes_fraction,
    exclude_genes = exclude_genes, exclude_gene_patterns = as.character(exclude_gene_patterns),
    seed = as.integer(seed)
  )
  return(res)
}

#' Compute forbidden genes
#' @param clean anndata object
#' @param suspect_gene_names gene names that might not be of value
#' @param suspect_gene_patterns gene patterns that might not be of value
#' @param seed seed
#'
#' @export
compute_forbidden_genes <- function(clean,
                                    suspect_gene_names = "",
                                    suspect_gene_patterns = "",
                                    seed = 123456) {
  reticulate::source_python(system.file("python", "metacells.py", package = "spacedeconv"))

  res <- compute_forbidden_genes(
    clean = clean,
    suspect_gene_names = suspect_gene_names,
    suspect_gene_patterns = suspect_gene_patterns,
    seed = as.integer(seed)
  )

  return(res)
}

#' extract forbidden genes from gene modules
#'
#' @param clean anndata object
#' @param forbidden_modules gene modules to remove
#'
#' @export
extract_forbidden_from_modules <- function(clean, forbidden_modules) {
  reticulate::source_python(system.file("python", "metacells.py", package = "spacedeconv"))

  res <- extract_forbidden_from_modules(clean = clean, forbidden_modules = forbidden_modules)

  return(res)
}

#' Compute Metacells
#' @param clean anndata object
#' @param forbidden_gene_names list of genes names that are not used as metacell base
#' @param cell_type_col cell type column of cleaned anndata, used for reannotation
#' @param abundance_score metacell celltype purity score
#' @export
compute_metacells <- function(clean, forbidden_gene_names, cell_type_col, abundance_score = 0.9) {
  reticulate::source_python(system.file("python", "metacells.py", package = "spacedeconv"))

  if (is.null(cell_type_col)) {
    stop("Please provide a cell type column name")
  }

  res <- compute_metacells(clean = clean, forbidden_gene_names = forbidden_gene_names)

  # reannotation
  internal <- res[[1]]
  metacell <- res[[2]]

  message("reannotating metacell types ...")

  celllist <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(celllist) <- c("metacell", "mostAbundant", "percentage")

  pb <- progress::progress_bar$new(total = length(rownames(metacell)))
  pb$tick(0)
  intDF <- internal$obs

  for (cell in rownames(metacell)) {
    # get the single cells added to this metacell
    tmp <- intDF[intDF$metacell == cell, ][, cell_type_col]
    # tmpcell <- list(list(tmp[cell_type_col]))
    tmpcell <- sort(table(unlist(tmp)), decreasing = T)[1]
    mostAbundant <- names(tmpcell)
    nMostAbundant <- unname(tmpcell)
    numberOfCells <- length(unlist(tmp))
    percent <- nMostAbundant / numberOfCells

    celllist[nrow(celllist) + 1, ] <- c(cell, mostAbundant, percent)
    pb$tick()
  }

  # print (metacell)

  metacell <- anndata_to_singlecellexperiment(metacell)
  SummarizedExperiment::assayNames(metacell) <- "metacell_counts" # renaming assay
  SummarizedExperiment::colData(metacell)$celltype <- celllist$mostAbundant # reannotation
  SummarizedExperiment::colData(metacell)$percentage <- as.numeric(celllist$percentage) # optional

  # filter for abundance score
  above90 <- celllist[celllist$percent >= abundance_score, ]$metacell
  message("Removing ", nrow(celllist) - length(above90), " metacell with abundance score under ", abundance_score)
  metacell <- metacell[, colnames(metacell) %in% above90]

  # scale
  tmp <- assay(metacell, "metacell_counts")
  for (i in 1:ncol(metacell)){
    tmp[, i] <- tmp[, i]/metacell$grouped[i]
  }

  assay(metacell, "metacell_scaled") <- tmp
  assay(metacell, "metacell_round") <- round(tmp)

  return(metacell)
}
