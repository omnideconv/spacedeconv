data("single_cell_data_2")
data("spatial_data_2")


test_that("SPOTlight deconvolution works", {
  # rest is already tested elsewhere

  # test that missing model is catched
  expect_error(object = SpaceDeconv::deconvolute(
    single_cell_obj = single_cell_data_2,
    signature = NULL,
    cell_type_col = "celltype_major",
    method = "spotlight",
    spatial_obj = spatial_data_2
  ), regexp = "Model is missing or null")

  # test that missing spatial_obj is catched
  expect_error(object = SpaceDeconv::deconvolute(
    single_cell_obj = single_cell_data_2,
    signature = NULL,
    cell_type_col = "celltype_major",
    method = "spotlight",
    spatial_obj = spatial_data_2
  ), regexp = "Model is missing or null")
})
