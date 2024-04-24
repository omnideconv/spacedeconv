deconvD <- readRDS(system.file("extdata", "deconv_dwls.rds", package = "spacedeconv"))
deconvQ <- readRDS(system.file("extdata", "deconv_quantiseq.rds", package = "spacedeconv"))

test_that("scatter plot works", {
  p <- plot_scatter(deconvD, "dwls_B.cells", deconvQ, "quantiseq_B.cell")

  expect_s3_class(p, "ggplot")
})

signature <- readRDS(system.file("extdata", "signature_dwls.rds", package = "spacedeconv"))

test_that("signature comparison works", {
  p <- compare_signatures(signature, signature)

  expect_s3_class(p, "ggplot")
})
