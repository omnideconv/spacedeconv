test_that("every registered deconvolution method has a real test case", {
  expect_setequal(names(method_cases), unname(deconvolution_methods))
})

# SPOTlight, Rectangle and DOT have dedicated files so they can be tested in isolation.
for (method in setdiff(names(method_cases), c("spotlight", "rectangle", "dot"))) {
  test_that(paste("real deconvolution:", method), {
    local_test_state()
    withr::local_dir(withr::local_tempdir())
    local_test_device()
    case <- method_cases[[method]]
    if (!is.null(case$package)) skip_if_not_installed(case$package)
    if (isTRUE(case$licensed)) local_cibersort()
    if (isTRUE(case$reference)) {
      fixture <- reference_mixture()
      spe <- fixture$spe
      sce <- fixture$sce
      common <- list(single_cell_obj = sce, spatial_obj = spe, method = method,
        cell_type_col = "celltype", batch_id_col = "sample_id")
      build_args <- list()
      deconv_args <- case$args
      if (method == "cell2location") {
        skip_if_not(reticulate::py_module_available("cell2location"), "Python cell2location is not installed")
        local_python_state()
        # Actual CPU training, deliberately small but not a one-epoch stub.
        build_args <- list(epochs = 100, gpu = FALSE, posterior_samples = 30)
        deconv_args <- list(epochs = 200, gpu = FALSE, posterior_samples = 30)
      }
      model <- do.call(build_model, c(common, build_args))
      if (method %in% c("spotlight", "spatialdwls", "cell2location")) expect_false(is.null(model))
      result <- do.call(deconvolute, c(common, list(signature = model, return_object = FALSE), deconv_args))
      expect_method_result(result, spe, method)
      fractions <- as.matrix(result)[rownames(fixture$truth), paste0(method, "_", colnames(fixture$truth)), drop = FALSE]
      # Strong, intentionally easy mixtures: most spots should recover the dominant type.
      expect_gte(mean(max.col(fractions, ties.method = "first") == max.col(fixture$truth)), .8)
    } else {
      spe <- if (isTRUE(case$mouse)) mouse_method_input() else human_method_input()
      result <- do.call(deconvolute, c(list(spatial_obj = spe, method = method, return_object = FALSE), case$args))
      expect_method_result(result, spe, method)
    }
  })
}

test_that("complete user workflow preserves data through plotting", {
  local_test_state()
  spe <- human_method_input()
  out <- build_and_deconvolute(small_sce(), spe, method = "estimate", cell_type_col = "celltype")
  expect_preserved(out, spe)
  columns <- available_results(out, "estimate")
  expect_gte(length(columns), 2)
  out <- aggregate_results(out, columns[1:2], name = "combined")
  expect_equal(out$combined, out[[columns[1]]] + out[[columns[2]]])
  expect_renderable(plot_spatial(out, "combined", density = FALSE))
})
