spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))
sce <- readRDS(system.file("testdata", "sce.rds", package = "spacedeconv"))


test_that("SPOTlight deconvolution works", {
  # test that missing model is catched
  expect_error(object = spacedeconv::deconvolute(
    single_cell_obj = sce,
    signature = NULL,
    cell_type_col = "celltype_major",
    method = "spotlight",
    spatial_obj = spe
  ), regexp = "Model is missing or null")

  # test that missing spatial_obj is catched
  expect_error(object = spacedeconv::deconvolute(
    single_cell_obj = sce,
    signature = NULL,
    cell_type_col = "celltype_major",
    method = "spotlight",
    spatial_obj = spe
  ), regexp = "Model is missing or null")
})
