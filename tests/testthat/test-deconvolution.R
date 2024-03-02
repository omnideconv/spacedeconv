spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))
sce <- readRDS(system.file("testdata", "sce.rds", package = "spacedeconv"))

testthat::test_that("QuanTIseq deconvolution works", {
  deconv <- deconvolute(spe, method = "quantiseq")

  # check that a column with "quantiseq" as token exists
  expect_true(any(grepl("quantiseq", colnames(colData(deconv)))))
})


test_that("Build_and_deconvolute works", {
  deconv <- build_and_deconvolute(spatial_obj = spe, method = "estimate", single_cell_obj = sce, cell_type_col = "celltype_major")

  expect_s4_class(deconv, "SpatialExperiment")
})
