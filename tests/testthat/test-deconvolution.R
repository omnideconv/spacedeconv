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

test_that("CARD deconvolution works", {
  deconv <- SpaceDeconv::deconvolute(spatial_obj = spatial_data_2,
                                     single_cell_obj = single_cell_data_2,
                                     signature = NULL,
                                     method="card",
                                     cell_type_col = "celltype_major",
                                     batch_id_col = "orig.ident")
  # test that deconvolution results are available
  expect_true(object = any(grepl(colnames(colData(deconv)), pattern = "card_")))
})

test_that("Mouse deconvolution works", {
  spatial_data_mouse <- spatial_data_2
  # mimic mouse data
  rownames(spatial_data_mouse) <- tools::toTitleCase(tolower(rownames(spatial_data_mouse)))

  deconv <- SpaceDeconv::deconvolute(spatial_data_mouse, method="mmcp_counter")

  expect_true(object = any(grepl(colnames(colData(deconv)), pattern = "mmcp_counter_")))

})
