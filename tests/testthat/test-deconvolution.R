spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))

testthat::test_that("QuanTIseq deconvolution works", {
  deconv = deconvolute(spe, method = "quantiseq")

  # check that a column with "quantiseq" as token exists
  expect_true(any(grepl("quantiseq", colnames(colData(deconv)))))
})
