spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))
sce <- readRDS(system.file("testdata", "sce.rds", package = "spacedeconv"))

spe <- spacedeconv::deconvolute(spe, method = "estimate")

test_that("plot_celltype executes without error", {
  expect_s3_class(plot_celltype(spe, cell_type = "estimate_tumor.purity"), "ggplot")
})

test_that("plot_umi_count executes without error", {
  expect_silent(plot_umi_count(spe))
})

test_that("plot_most_abundant executes without error", {
  expect_silent(plot_most_abundant(spe, method = "estimate"))
})

test_that("plot_gene executes without error", {
  expect_silent(plot_gene(spe, gene = "MOV10"))
})
