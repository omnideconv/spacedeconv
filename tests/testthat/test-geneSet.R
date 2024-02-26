spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))


test_that("gene_set_score stops with invalid spe input", {
  expect_error(gene_set_score(NULL, genes = c("CD58", "DR1")), "spe must be a SpatialExperiment object.")
})

test_that("gene_set_score stops with invalid genes input", {
  expect_error(gene_set_score(spe, genes = NULL), "genes must be a non-null vector of gene names.")
  expect_error(gene_set_score(spe, genes = list("CD58", "DR1")), "genes must be a non-null vector of gene names.")
})

test_that("gene_set_score stops with invalid name input", {
  expect_error(gene_set_score(spe, genes = c("CD58", "DR1"), name = NULL), "name must be a non-empty single string.")
})

test_that("gene_set_score calculates score with default assay", {
  modified_spe <- gene_set_score(spe, genes = c("CD58", "DR1"), name = "testScore", assay = "counts")
  expect_true("testScore" %in% colnames(colData(modified_spe)))
})
