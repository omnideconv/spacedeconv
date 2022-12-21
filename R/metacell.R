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
  # reticulate::source_python("./inst/python/metacells.py")
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
#' @export
compute_metacells <- function(clean, forbidden_gene_names) {
  reticulate::source_python(system.file("python", "metacells.py", package = "spacedeconv"))

  res <- compute_metacells(clean = clean, forbidden_gene_names = forbidden_gene_names)
  return(res)
}
