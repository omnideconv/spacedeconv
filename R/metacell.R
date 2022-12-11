# library(reticulate)
# library(SingleCellExperiment)
#
# mc = reticulate::import("metacells")
# ad = reticulate::import("anndata")
#
# seurat = readRDS("~/subsetData/breast_cancer_seurat_object.rds")
#
# sce = Seurat::as.SingleCellExperiment(seurat)
#
#
# # observeration = cell
# # variables = gene
# #anndata = ad$AnnData(X = t(SingleCellExperiment::counts(sce)),
#   #                   var = as.data.frame(rowData(sce)), obs = as.data.frame(colData(sce)))
#
#
#
# # variables
# excluded_gene_names = list('IGHMBP2', 'IGLL1', 'IGLL5', 'IGLON5', 'NEAT1', 'TMSB10', 'TMSB4X')
# excluded_gene_patterns = list('MT-.*')
#
#
#
# mc$pl$analyze_clean_genes(anndata, excluded_gene_names = excluded_gene_names,
#                           excluded_gene_patterns = excluded_gene_patterns,
#                           random_seed = 123456)
#
# anndata

#' Clean Genes and Cells with Metacell2
#'
#' @param anndata anndata object containing Single Cell Expression Data
#' @param properly_sampled_min_cell_total testtesttest
#' @param properly_sampled_max_cell_total testtesttest
#' @param properly_sampled_max_excluded_genes_fraction testtesttest
#' @param exclude_genes test
#' @param exclude_gene_patterns test
#' @param seed seed
#'
#' @export
clean_genes_and_cells <- function(anndata, properly_sampled_min_cell_total = 800,
                                  properly_sampled_max_cell_total = 8000,
                                  properly_sampled_max_excluded_genes_fraction = 0.1,
                                  exclude_genes = "", exclude_gene_patterns = "", seed = 123456){

}

#' Compute forbidden genes
#' @param clean anndata object
#' @param suspect_gene_names test
#' @param suspect_gene_patterns test
#' @param seed seed
#'
#' @export
compute_forbidden_genes <- function (clean, suspect_gene_names="", suspect_gene_patterns="", seed = 123456){

}

#' extract forbidden genes from gene modules
#'
#' @param clean anndata object
#' @param forbidden_modules test
#'
#' @export
extract_forbidden_from_modules <- function (clean, forbidden_modules){

}

#' Compute Metacells
#' @param clean anndata object
#' @param forbidden_gene_names test
compute_metacells <- function (clean, forbidden_gene_names){

}
