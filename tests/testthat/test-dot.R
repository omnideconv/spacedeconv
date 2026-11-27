test_that("DOT runs through the public workflow without attaching DOTr", {
  local_test_state()
  withr::local_dir(withr::local_tempdir())
  local_test_device()
  expect_true(requireNamespace("DOTr", quietly = TRUE))
  expect_false("package:DOTr" %in% search())

  fixture <- reference_mixture()
  expect_null(build_model(fixture$sce, method = "dot", cell_type_col = "celltype"))
  result <- deconvolute(
    spatial_obj = fixture$spe,
    single_cell_obj = fixture$sce,
    cell_type_col = "celltype",
    method = "dot",
    return_object = FALSE
  )
  expect_method_result(result, fixture$spe, "dot")
  fractions <- as.matrix(result)[rownames(fixture$truth),
    paste0("dot_", colnames(fixture$truth)), drop = FALSE]
  expect_gte(mean(max.col(fractions, ties.method = "first") == max.col(fixture$truth)), .8)
  expect_false("package:DOTr" %in% search())
})
