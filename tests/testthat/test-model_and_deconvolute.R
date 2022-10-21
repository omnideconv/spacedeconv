data("single_cell_data_2")
data("spatial_data_2")

test_that("RCTD Model and deconvolution in one step works", {
  deconv <- spacedeconv::build_and_deconvolute(
    single_cell_obj = single_cell_data_2,
    spatial_obj = spatial_data_2,
    method = "rctd",
    cell_type_col = "celltype_major",
    return_object = FALSE
  )

  expect_equal(
    info = "deconvolution results for all cell types in single cell reference available",
    obj = sort(colnames(deconv)), expected = sort(paste0("rctd_", unique(single_cell_data_2$celltype_major)))
  )
})

test_that("SPOTlight model and deconvolution in one step works", {
  deconv <- spacedeconv::build_and_deconvolute(
    single_cell_obj = single_cell_data_2,
    spatial_obj = spatial_data_2,
    method = "spotlight",
    cell_type_col = "celltype_major",
    return_object = FALSE
  )

  expect_equal(
    info = "deconvolution results for all cell types in single cell reference available",
    obj = sort(colnames(deconv)), expected = sort(paste0("spotlight_", unique(single_cell_data_2$celltype_major)))
  )
})

# test_that("CARD model and deconvolution in one step works", {
#   deconv <- spacedeconv::build_and_deconvolute(
#     single_cell_obj = single_cell_data_2,
#     spatial_obj = spatial_data_2,
#     method = "card",
#     cell_type_col = "celltype_major",
#     return_object = FALSE
#   )
#
#   expect_equal(
#     info = "deconvolution results for all cell types in single cell reference available",
#     obj = sort(colnames(deconv)), expected = sort(unique(single_cell_data_2$celltype_major))
#   )
# })


# currently the only omnideconv method
# test_that("DWLS model and deconvolution in one step works", {
#   deconv <- spacedeconv::build_and_deconvolute(
#     single_cell_obj = single_cell_data_2,
#     spatial_obj = spatial_data_2,
#     method = "dwls",
#     cell_type_col = "celltype_major",
#     return_object = FALSE,
#     dwls_method = "mast_optimized"
#   )
#
#   expect_equal(
#     info = "deconvolution results for all cell types in single cell reference available",
#     obj = sort(colnames(deconv)), expected = sort(unique(single_cell_data_2$celltype_major))
#   )
# })

# currently the only immunedeconv method
test_that("quanTIseq model and deconvolution in one step works", {
  deconv <- spacedeconv::build_and_deconvolute(
    single_cell_obj = single_cell_data_2,
    spatial_obj = spatial_data_2,
    method = "quantiseq",
    cell_type_col = "celltype_major",
    return_object = FALSE,
    dwls_method = "mast_optimized"
  )

  expect_gt(object = dim(deconv)[1], expected = 1)

  expect_gt(object = dim(deconv)[2], expected = 1)
})
