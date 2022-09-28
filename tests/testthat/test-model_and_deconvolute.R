data("single_cell_data_2")
data("spatial_data_2")

test_that("RCTD Model and deconvolution in one step works", {
  deconv <- SpaceDeconv::build_and_deconvolute(
    single_cell_obj = single_cell_data_2,
    spatial_obj = spatial_data_2,
    method = "rctd",
    cell_type_col="celltype_major"
  )



})
