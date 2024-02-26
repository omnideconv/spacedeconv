spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))
sce <- readRDS(system.file("testdata", "sce.rds", package = "spacedeconv"))

spe <- spacedeconv::deconvolute(spe, method = "epic")

test_that("plot_celltype executes without error", {
  expect_silent(plot_celltype(spe, cell_type = "epic_NK.cell"))
})

test_that("plot_umi_count executes without error", {
  expect_silent(plot_umi_count(spe))
})

test_that("plot_most_abundant executes without error", {
  expect_silent(plot_most_abundant(spe, method = "epic"))
})

test_that("plot_comparison executes without error", {
  expect_silent(plot_comparison(spe, cell_type_1 = "epic_B.cell", cell_type_2 = "epic_NK.cell"))
})

test_that("plot_gene executes without error", {
  expect_silent(plot_gene(spe, gene = "MOV10"))
})

