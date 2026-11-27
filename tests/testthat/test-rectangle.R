test_that("Rectangle validates expression inputs before calling Python", {
  sce <- small_sce()
  spe <- small_spe()
  expect_error(rectangle_prepare_bulks(spe, "missing"), "assay is not available")
  expect_error(rectangle_prepare_single_cell(sce, "celltype", "missing"), "assay is not available")
  expect_error(build_model_rectangle(sce, cell_type_col = "missing"), "can't be found")
  expect_error(deconvolute_rectangle(spe), "signature")
  expect_error(rectangle_check_cpus(1.5), "positive integer")
  expect_error(rectangle_check_cpus(0), "positive integer")
  assay(sce, "counts")[1, 1] <- .5
  expect_error(rectangle_prepare_single_cell(sce, "celltype", "counts"), "raw integer counts")
  assay(spe, "counts")[1, 1] <- NA_real_
  expect_error(rectangle_prepare_bulks(spe, "counts"), "finite, non-negative")
  spe <- small_spe()
  assay(spe, "counts")[, 1] <- 0
  expect_error(rectangle_prepare_bulks(spe, "counts"), "positive library sizes")
  spe <- small_spe()
  colnames(spe)[2] <- colnames(spe)[1]
  expect_error(rectangle_prepare_bulks(spe, "counts"), "unique, non-empty")
})

test_that("Rectangle aligns results by spot ID and rejects incompatible results", {
  spe <- small_spe()
  expected <- matrix(rep(c(.2, .8), each = ncol(spe)), ncol(spe),
    dimnames = list(colnames(spe), c("A", "B")))
  expected[, 1] <- seq(.1, .6, by = .1)
  result <- expected[rev(seq_len(nrow(expected))), , drop = FALSE]
  testthat::local_mocked_bindings(rectangle_python = function() {
    list(py_deconvolute_rectangle = function(...) result)
  })
  actual <- deconvolute_rectangle(spe, signature = list())
  colnames(expected) <- paste0("rectangle_", colnames(expected))
  expect_equal(actual, expected)
  rownames(result)[1] <- "wrong_spot"
  expect_error(deconvolute_rectangle(spe, signature = list()), "spot IDs")
})

test_that("real Rectangle workflow recovers simple mixtures and preserves the spatial object", {
  skip_if_not(reticulate::py_module_available("rectanglepy"), "Python rectanglepy is not installed")
  local_test_state()
  local_test_device()
  fixture <- reference_mixture()
  signature <- build_model(fixture$sce, spatial_obj = fixture$spe,
    method = "rectangle", cell_type_col = "celltype", n_cpus = 1L,
    optimize_cutoffs = FALSE)
  expect_true(reticulate::is_py_object(signature))
  result <- deconvolute(fixture$spe, signature = signature, method = "rectangle",
    n_cpus = 1L, return_object = FALSE)
  expect_method_result(result, fixture$spe, "rectangle")
  expect_true(all(result >= -1e-8))
  fractions <- as.matrix(result)[rownames(fixture$truth), paste0("rectangle_", colnames(fixture$truth)), drop = FALSE]
  expect_gte(mean(max.col(fractions, ties.method = "first") == max.col(fixture$truth)), .8)
  out <- deconvolute(fixture$spe, signature = signature, method = "Rectangle", n_cpus = 1L)
  expect_preserved(out, fixture$spe)
  expect_equal(as.matrix(colData(out)[, colnames(result), drop = FALSE]), result, ignore_attr = TRUE)
  expect_false(exists("py_build_rectangle_signatures", envir = .GlobalEnv, inherits = FALSE))
})
