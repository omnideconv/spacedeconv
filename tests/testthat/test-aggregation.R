deconvD <- readRDS(system.file("extdata", "deconv_dwls.rds", package = "spacedeconv"))

test_that("aggregation works", {
  agg <- aggregate_results(deconvD, cell_types = c("dwls_Normal.Epithelial", "dwls_Cancer.Epithelial"), name = "test")

  expect_contains(names(colData(agg)), "test")
})
