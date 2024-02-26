spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))
sce <- readRDS(system.file("testdata", "sce.rds", package = "spacedeconv"))

testthat::test_that("QuanTIseq deconvolution works", {
  deconv <- deconvolute(spe, method = "quantiseq")

  # check that a column with "quantiseq" as token exists
  expect_true(any(grepl("quantiseq", colnames(colData(deconv)))))
})



testthat::test_that("SCDC deconvolution works", {
  deconv = deconvolute(spe, method = "scdc", single_cell_obj = sce, cell_type_col = "celltype_major", batch_id_col = "orig.ident")
  expect_true(any(grepl("scdc", available_results(deconv))))
})
