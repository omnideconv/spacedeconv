test_that("method registry is unambiguous and signature-free methods return NULL", {
  expect_false(anyDuplicated(unname(deconvolution_methods)) > 0)
  expect_true(all(nzchar(names(deconvolution_methods))))
  for (method in c(unname(first_gen), "rctd", "card", "dot")) {
    expect_null(build_model(small_sce(), cell_type_col = "celltype", method = method))
  }
  expect_error(build_model(NULL, method = "rctd"), "missing or null")
  expect_error(build_model(small_sce(), method = "unknown"), "not supported")
  expect_error(deconvolute(NULL), "missing or null")
  expect_error(deconvolute(small_spe(), method = "unknown"), "not recognized")
})

test_that("documented display names work for model building", {
  expect_null(build_model(small_sce(), method = "RCTD", cell_type_col = "celltype"))
})

test_that("the convenience wrapper forwards the model, assays and return mode", {
  model <- matrix(1, 2, 2)
  raw <- matrix(.5, 6, 2)
  received <- NULL
  local_mocked_bindings(
    build_model = function(...) model,
    deconvolute = function(...) { received <<- list(...); raw })
  out <- build_and_deconvolute(small_sce(), small_spe(), method = "rctd",
    cell_type_col = "celltype", assay_sc = "reference_counts",
    assay_sp = "spatial_counts", return_object = FALSE, n_cores = 1)
  expect_identical(out, raw)
  expect_identical(received$signature, model)
  expect_identical(received$assay_sc, "reference_counts")
  expect_identical(received$assay_sp, "spatial_counts")
  expect_false(received$return_object)
  expect_equal(received$n_cores, 1)
})
