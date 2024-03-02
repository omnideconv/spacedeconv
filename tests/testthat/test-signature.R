spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))
sce <- readRDS(system.file("testdata", "sce.rds", package = "spacedeconv"))



test_that("Immunedeconv models works", {
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "mcp_counter"), info = "Signature is null (which it should be)")
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "epic"), info = "Signature is null (which it should be)")
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "quantiseq"), info = "Signature is null (which it should be)")
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "xcell"), info = "Signature is null (which it should be)")
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "cibersort"), info = "Signature is null (which it should be)")
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "cibersort_abs"), info = "Signature is null (which it should be)")
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "timer"), info = "Signature is null (which it should be)")
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "consensus_tme"), info = "Signature is null (which it should be)")
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "abis"), info = "Signature is null (which it should be)")
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "estimate"), info = "Signature is null (which it should be)")
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "mmcp_counter"), info = "Signature is null (which it should be)")
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "seqimmucc"), info = "Signature is null (which it should be)")
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "dcq"), info = "Signature is null (which it should be)")
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "base"), info = "Signature is null (which it should be)", )
})



# Test for input validation in build_model
test_that("build_model requires non-null single_cell_obj", {
  expect_error(build_model(NULL, method = "rctd"), "Parameter 'single_cell_obj' missing or null, but is required")
})

test_that("build_model requires supported method", {
  expect_error(build_model(sce, method = "unsupported_method"), "Parameter 'method' is null or not supported")
})

# Test for functionality of deconvolute
test_that("deconvolute requires non-null spatial_obj", {
  expect_error(deconvolute(NULL), "Parameter 'spatial_obj' is missing or null, but is required.")
})

test_that("Second-gen signature creation with SCDC works", {
  signature <- build_model(sce, cell_type_col = "celltype_major", method = "scdc", verbose = T, batch_id_col = "orig.ident")
  expect_equal(ncol(signature), expected = 4)
  expect_type(signature, "double")
})


test_that("Second-gen deconvolution with SCDC works", {
  deconv <- deconvolute(spe, signature, method = "scdc", single_cell_obj = sce, cell_type_col = "celltype_major", batch_id_col = "orig.ident")
  expect_true("scdc_Endothelial" %in% colnames(colData(deconv)))
  expect_s4_class(deconv, "SpatialExperiment")
})
