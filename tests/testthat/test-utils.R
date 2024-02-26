sce <- readRDS(system.file("testdata", "sce.rds", package = "spacedeconv"))
library(SingleCellExperiment)

# Test the subsetSCE function
testthat::test_that("subsetSCE subsets the SingleCellExperiment correctly", {
  # Call the subsetSCE function
  result <- subsetSCE(sce, cell_type_col = "celltype_major", scenario = "even", ncells = 500, notEnough = "all", seed = 12345)

  # Check if the resulting object is a SingleCellExperiment
  expect_true(inherits(result, "SingleCellExperiment"))

  # Check the number of cells in the resulting objectpak:
  expect_equal(ncol(result), 250)
})



testthat::test_that("Normalization works", {
  normalized <- spacedeconv::normalize(sce, method = "cpm", assay = "counts")

  expect_true("cpm" %in% assayNames(normalized))
})



