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
  cli::cli_rule(left = "metacell")

  cli::cli_progress_step("Cleaning genes and cells", msg_done = "Cleaned genes and cells")

  #init_python()
  #metacells_checkload()

  reticulate::source_python(system.file("python", "metacells.py", package = "spacedeconv"))

  res <- clean_genes_and_cells(anndata,
    properly_sampled_min_cell_total = as.integer(properly_sampled_min_cell_total),
    properly_sampled_max_cell_total = as.integer(properly_sampled_max_cell_total),
    properly_sampled_max_excluded_genes_fraction = properly_sampled_max_excluded_genes_fraction,
    exclude_genes = exclude_genes, exclude_gene_patterns = as.character(exclude_gene_patterns),
    seed = as.integer(seed)
  )

  cli::cli_progress_done()

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
  cli::cli_rule(left = "metacell")

  cli::cli_progress_step("Computing forbidden genes", msg_done = "Computed forbidden genes")

  #init_python()
  #metacells_checkload()

  reticulate::source_python(system.file("python", "metacells.py", package = "spacedeconv"))

  res <- compute_forbidden_genes(
    clean = clean,
    suspect_gene_names = suspect_gene_names,
    suspect_gene_patterns = suspect_gene_patterns,
    seed = as.integer(seed)
  )

  cli::cli_progress_done()

  return(res)
}

#' extract forbidden genes from gene modules
#'
#' @param clean anndata object
#' @param forbidden_modules gene modules to remove
#'
#' @export
extract_forbidden_from_modules <- function(clean, forbidden_modules) {
  cli::cli_rule(left = "metacell")

  cli::cli_progress_step("Extracting forbidden genes", msg_done = "Extracted forbidden genes")

  #init_python()
  #metacells_checkload()

  reticulate::source_python(system.file("python", "metacells.py", package = "spacedeconv"))

  res <- extract_forbidden_from_modules(clean = clean, forbidden_modules = forbidden_modules)

  cli::cli_progress_done()

  return(res)
}

#' Compute Metacells
#' @param clean anndata object
#' @param forbidden_gene_names list of genes names that are not used as metacell base
#' @param cell_type_col cell type column of cleaned anndata, used for reannotation
#' @param target_size target UMI count per metacells
#' @param abundance_score metacell celltype purity score
#' @param n_cores number of cores to use
#' @param seed seed
#' @export
compute_metacells <- function(clean, forbidden_gene_names, cell_type_col, target_size = 160000, abundance_score = 0.9, n_cores = NULL, seed = 12345) {
  cli::cli_rule(left = "metacell")

  cli::cli_progress_step("Computing metacells", msg_done = "computed metacells")

  #init_python()
  #metacells_checkload()

  reticulate::source_python(system.file("python", "metacells.py", package = "spacedeconv"))

  if (is.null(cell_type_col)) {
    stop("Please provide a cell type column name")
  }

  if (is.numeric(n_cores)) {
    Sys.setenv(METACELLS_PROCESSORS_COUNT = n_cores)
  }

  res <- compute_metacells(clean = clean, forbidden_gene_names = forbidden_gene_names, target_size = as.integer(target_size), seed = as.integer(seed))

  # reannotation
  internal <- res[[1]]
  metacell <- res[[2]]

  # cli::cli_progress_step("Reannotating cell types", msg_done = "Reannotated cell types")

  celllist <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(celllist) <- c("metacell", "mostAbundant", "percentage")

  # pb <- progress::progress_bar$new(total = length(rownames(metacell)))
  # pb$tick(0)
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
    # pb$tick()
  }


  # print (metacell)

  metacell <- anndata_to_singlecellexperiment(metacell)
  SummarizedExperiment::assayNames(metacell) <- "counts" # renaming assay
  SummarizedExperiment::colData(metacell)$celltype <- celllist$mostAbundant # reannotation
  SummarizedExperiment::colData(metacell)$percentage <- as.numeric(celllist$percentage) # optional

  # filter for abundance score
  above90 <- celllist[celllist$percent >= abundance_score, ]$metacell
  message("Removing ", nrow(celllist) - length(above90), " metacell with abundance score under ", abundance_score)
  metacell <- metacell[, colnames(metacell) %in% above90]

  # scale
  tmp <- assay(metacell, "counts")
  for (i in 1:ncol(metacell)) {
    tmp[, i] <- tmp[, i] / metacell$grouped[i]
  }

  assay(metacell, "scaled") <- tmp
  assay(metacell, "round") <- round(tmp)

  # make them sparse
  assay(metacell, "counts") <- Matrix::Matrix(assay(metacell, "counts"), sparse = TRUE)
  assay(metacell, "scaled") <- Matrix::Matrix(assay(metacell, "scaled"), sparse = TRUE)
  assay(metacell, "round") <- Matrix::Matrix(assay(metacell, "round"), sparse = TRUE)

  cli::cli_progress_done()

  return(metacell)
}
