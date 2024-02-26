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


test_that("CARD signature works", {
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "card"), info = "Signature is null (what it should be)")
})

test_that("Bisque signature works", {
  expect_null(object = spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "bisque"), info = "Signature is null (what it should be)")
})
