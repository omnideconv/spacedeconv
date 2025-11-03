spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))
sce <- readRDS(system.file("testdata", "sce.rds", package = "spacedeconv"))


test_that("convert human to mouse with valid gene symbols", {
  human_genes <- c("BRCA1", "TP53") # Assuming these are valid human gene symbols
  result <- convert_human_to_mouse(human_genes)
  expect_type(result, "list")
  expect_equal(nrow(result), length(human_genes))
  expect_true(all(result$Human_symbol %in% human_genes))
})


ad <- spe_to_ad(spe, assay = "counts")

test_that("spe_to_ad returns an AnnData object", {
  expect_true("AnnDataR6" %in% class(ad))
})


test_that("spe_to_ad handles NULL input", {
  expect_error(spe_to_ad(NULL))
})

test_that("spe_to_ad handles non-existent assay", {
  expect_error(spe_to_ad(spe, assay = "nonexistent"))
})
