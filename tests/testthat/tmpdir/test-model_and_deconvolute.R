spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))
sce <- readRDS(system.file("testdata", "sce.rds", package = "spacedeconv"))




test_that("quanTIseq model and deconvolution in one step works", {
  deconv <- spacedeconv::build_and_deconvolute(
    single_cell_obj = sce,
    spatial_obj = spe,
    method = "quantiseq",
    cell_type_col = "celltype_major",
    return_object = FALSE
  )

  expect_gt(object = dim(deconv)[1], expected = 1)

  expect_gt(object = dim(deconv)[2], expected = 1)
})
