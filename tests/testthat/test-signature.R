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



test_that("Bisque signature creation works", {
  signature <- spacedeconv::build_model(single_cell_obj = sce, cell_type_col = "celltype_major", method = "bisque")
  expect_null(object = signature, info = "Bisque signature is null (which it should be)")
})




test_that("MOMF signature creation works", {
  signature <- spacedeconv::build_model(
    single_cell_obj = sce,
    cell_type_col = "celltype_major",
    method = "momf",
    spatial_obj = spe
  )
  expect_equal(
    info = "Signature matrix has the same number of celltypes as unique celltypes as reference data",
    object = ncol(signature),
    expected = length(unique(sce$celltype_major))
  )

  expect(
    info = "signature matrix genes are contained in the single cell reference",
    ok = sum(rownames(signature) %in% rownames(sce)) == dim(signature)[[1]],
    failure_message = "Not all signature genes are included in the single cell reference data"
  )
})




test_that("SCDC signature creation works", {
  signature <- spacedeconv::build_model(
    single_cell_obj = sce,
    cell_type_col = "celltype_major",
    method = "scdc", batch_id_col = "orig.ident"
  )
  expect_equal(
    info = "Signature matrix contains the same amount of unique celltypes as the single cell reference",
    object = ncol(signature),
    expected = length(unique(sce$celltype_major))
  )
})

test_that("RCTD signature creation works", {
  signature <- spacedeconv::build_model(
    single_cell_obj = sce,
    cell_type_col = "celltype_major",
    method = "rctd"
  )
  expect_null(object = signature)
})


test_that("CARD signature creation works", {
  signature <- spacedeconv::build_model(
    single_cell_obj = sce,
    cell_type_col = "celltype_major",
    method = "card"
  )
  expect_null(
    info = "CARD signature is null (what it should be)",
    object = signature
  )
})
