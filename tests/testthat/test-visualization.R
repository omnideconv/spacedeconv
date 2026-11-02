spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))
spe <- deconvolute(spe, method = "estimate")


test_that("Plot celltype works", {
  p <- plot_celltype(spe, cell_type = "estimate_estimate.score", show_image = TRUE)

  expect_s3_class(p, "ggplot")
})


test_that("Plot Gene works", {
  p <- plot_gene(spe, gene = "CCN1", transform_scale = "log")

  expect_s3_class(p, "ggplot")
})


test_that("Plot Umi Count works", {
  p <- plot_umi_count(spe, smooth = TRUE)

  expect_s3_class(p, "ggplot")
})


test_that("Plot_post_abundant works", {
  p <- plot_most_abundant(spe, method = "estimate")

  expect_s3_class(p, "ggplot")
})
