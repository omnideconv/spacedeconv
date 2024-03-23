deconvD <- readRDS(system.file("extdata", "deconv_dwls.rds", package = "spacedeconv"))

test_that("spatialcorr works", {
  res <- spatialcorr(deconvD, "dwls")

  expect_type(res, "list")
})
