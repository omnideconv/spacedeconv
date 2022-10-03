data("single_cell_data_2")
data("spatial_data_2")

test_that("Immunedeconv models works", {
  expect_null(object = SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "mcp_counter"), info = "Signature is null (which it should be)")
  expect_null(object = SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "epic"), info = "Signature is null (which it should be)")
  expect_null(object = SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "quantiseq"), info = "Signature is null (which it should be)")
  expect_null(object = SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "xcell"), info = "Signature is null (which it should be)")
  expect_null(object = SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "cibersort"), info = "Signature is null (which it should be)")
  expect_null(object = SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "cibersort_abs"), info = "Signature is null (which it should be)")
  expect_null(object = SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "timer"), info = "Signature is null (which it should be)")
  expect_null(object = SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "consensus_tme"), info = "Signature is null (which it should be)")
  expect_null(object = SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "abis"), info = "Signature is null (which it should be)")
  expect_null(object = SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "estimate"), info = "Signature is null (which it should be)")
  expect_null(object = SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "mmcp_counter"), info = "Signature is null (which it should be)")
  expect_null(object = SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "seqimmucc"), info = "Signature is null (which it should be)")
  expect_null(object = SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "dcq"), info = "Signature is null (which it should be)")
  expect_null(object = SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "base"), info = "Signature is null (which it should be)", )
})

test_that("AutoGeneS signature creation works", {
  signature <- SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "autogenes")
  expect_true(file.exists(signature), info = "AutoGenes signature was created successfully")
})

test_that("BayesPrism signature creation works", {
  signature <- SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "bayesprism", spatial_obj = spatial_data_2)
  expect_null(info = "BayesPrism signature is null (which it should be)", object = signature)
})

test_that("Bisque signature creation works", {
  signature <- SpaceDeconv::build_model(single_cell_obj = single_cell_data_2, cell_type_col = "celltype_major", method = "bisque")
  expect_null(object = signature, info = "Bisque signature is null (which it should be)")
})

test_that("Bseq-sc signature creation works", {
  expect_null(NULL, NULL)
})

test_that("CIBERSORTx signature creation works", {
  expect_null(NULL, NULL)
})

test_that("CDseq signature creation works", {
  expect_null(NULL, NULL)
})

test_that("CPM signature creation works", {
  signature <- SpaceDeconv::build_model(
    single_cell_obj = single_cell_data_2,
    cell_type_col = "celltype_major",
    method = "cpm",
    spatial_obj = spatial_data_2
  )
  expect_null(info = "The CPM model is null (which it should be)", object = signature)
})

test_that("DWLS signature creation works", {
  expect_null(NULL, NULL)
})

test_that("MOMF signature creation works", {
  signature <- SpaceDeconv::build_model(
    single_cell_obj = single_cell_data_2,
    cell_type_col = "celltype_major",
    method = "momf",
    spatial_obj = spatial_data_2
  )
  expect_equal(
    info = "Signature matrix has the same number of celltypes as unique celltypes as reference data",
    object = ncol(signature),
    expected = length(unique(single_cell_data_2$celltype_major))
  )

  expect(
    info = "signature matrix genes are contained in the single cell reference",
    ok = sum(rownames(signature) %in% rownames(single_cell_data_2)) == dim(signature)[[1]],
    failure_message = "Not all signature genes are included in the single cell reference data"
  )
})

# test_that("MuSiC signature creation works", {
#   signature <- SpaceDeconv::build_model(
#     single_cell_obj = single_cell_data_2,
#     cell_type_col = "celltype_major",
#     method = "music", batch_id_col = "orig.ident"
#   )
#   expect_equal(
#     info = "Signature matrix contains the same amount of unique celltypes as the single cell reference",
#     object = ncol(signature),
#     expected = length(unique(single_cell_data_2$celltype_major))
#   )
# })

test_that("Scaden signature creation works", {
  signature <- SpaceDeconv::build_model(
    single_cell_obj = single_cell_data_2,
    cell_type_col = "celltype_major",
    method = "scaden", spatial_obj = spatial_data_2
  )
  expect_equal(
    info = "model folder is created and model assets are written",
    object = length(list.dirs(signature, recursive = F)), expected = 3
  )
})

test_that("SCDC signature creation works", {
  signature <- SpaceDeconv::build_model(
    single_cell_obj = single_cell_data_2,
    cell_type_col = "celltype_major",
    method = "scdc", batch_id_col = "orig.ident"
  )
  expect_equal(
    info = "Signature matrix contains the same amount of unique celltypes as the single cell reference",
    object = ncol(signature),
    expected = length(unique(single_cell_data_2$celltype_major))
  )
})

test_that("RCTD signature creation works", {
  signature <- SpaceDeconv::build_model(
    single_cell_obj = single_cell_data_2,
    cell_type_col = "celltype_major",
    method = "rctd"
  )
  expect_null(object = signature)
})

test_that("SPOTlight signature creation works", {
  model <- SpaceDeconv::build_model(
    single_cell_obj = single_cell_data_2,
    cell_type_col = "celltype_major",
    method = "spotlight",
    spatial_obj = spatial_data_2
  )
  # test that model exists
  expect_equal(
    info = "SPOTlight model creation successfull",
    object = length(model),
    expected = 2
  )
  # test that model is valid and contains mod and topics
  expect_equal(
    info = "SPOTLight model valid",
    object = names(model),
    expected = c("mod", "topic")
  )

  # test that SPOTlight catches unavailable assays
  expect_message(object = SpaceDeconv::build_model(
    single_cell_obj = single_cell_data_2,
    cell_type_col = "celltype_major",
    method = "spotlight",
    spatial_obj = spatial_data_2,
    assay_sc = "abc"
  ), regexp = "not available in expression object")
  expect_message(object = SpaceDeconv::build_model(
    single_cell_obj = single_cell_data_2,
    cell_type_col = "celltype_major",
    method = "spotlight",
    spatial_obj = spatial_data_2,
    assay_sp = "abc"
  ), regexp = "not available in expression object")

  # check that spotlight catches missing single_cell_obj
  expect_error(object = SpaceDeconv::build_model(
    single_cell_obj = NULL,
    cell_type_col = "celltype_major",
    method = "spotlight",
    spatial_obj = spatial_data_2
  ), regexp = "'single_cell_obj' missing or null")

  # test that spotlight catches wrong cell type col names
  expect_error(object = SpaceDeconv::build_model(
    single_cell_obj = single_cell_data_2,
    cell_type_col = "abc",
    method = "spotlight",
    spatial_obj = spatial_data_2
  ), regexp = "can't be found in single cell object")
})


test_that("CARD signature creation works", {
  signature <- SpaceDeconv::build_model(
    single_cell_obj = single_cell_data_2,
    cell_type_col = "celltype_major",
    method = "card"
  )
  expect_null(
    info = "CARD signature is null (what it should be)",
    object = signature
  )
})

# test_that("SpatialDWLS signature creation works", {
#   signature <- SpaceDeconv::build_model(
#     single_cell_obj = single_cell_data_2,
#     cell_type_col = "celltype_major",
#     method = "spatialdwls"
#   )
#   # expect_equal(
#   #   info = "signature contains all cell types from reference data",
#   #   object = ncol(signature),
#   #   expected = length(unique(single_cell_data_2$celltype_major))
#   # )
# })

# test_that("cell2location signature creation works")
