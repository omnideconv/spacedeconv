test_that("SPOTlight requires explicit markers", {
  expect_error(build_model(small_sce(), spatial_obj = small_spe(),
    cell_type_col = "celltype", method = "spotlight"), "SPOTlight requires 'markers'")
  expect_error(build_and_deconvolute(small_sce(), small_spe(),
    cell_type_col = "celltype", method = "spotlight"), "SPOTlight requires 'markers'")
  expect_error(deconvolute(small_spe(), method = "spotlight"), "Model is missing")
  expect_error(build_model(small_sce(), cell_type_col = "celltype",
    method = "spotlight"), "spatial")
})

test_that("SPOTlight rejects invalid markers before training", {
  call <- function(markers) build_model_spotlight(small_sce(), "celltype", small_spe(), markers = markers)
  good <- data.frame(gene = c("G1", "G2"), cluster = c("A", "B"), mean.AUC = c(.9, .95))
  for (bad in list(list(), good[FALSE, ], good[, 1:2])) {
    expect_error(call(bad), "non-empty data frame")
  }
  for (weight in list(NA_real_, Inf, -1, 0, 1.1, "0.9")) {
    bad <- good; bad$mean.AUC <- weight
    expect_error(call(bad), "finite numeric values")
  }
  for (label in c(NA_character_, "", " ")) {
    bad <- good; bad$gene[1] <- label
    expect_error(call(bad), "non-empty strings")
  }
  expect_error(call(rbind(good, good[1, ])), "duplicate")
  expect_error(call(good[1, , drop = FALSE]), "match all cell types")
  bad <- good; bad$cluster[1] <- "unknown"
  expect_error(call(bad), "match all cell types")
  bad <- good; bad$gene[1] <- "absent"
  expect_error(call(bad), "present in both")
})

test_that("SPOTlight forwards supplied markers and training options", {
  received <- NULL
  local_mocked_bindings(trainNMF = function(...) { received <<- list(...); "model" },
    .package = "SPOTlight")
  markers <- S4Vectors::DataFrame(gene = c("G1", "G2"),
    cluster = factor(c("A", "B")), mean.AUC = c(.9, .95))
  out <- build_model(small_sce(), spatial_obj = small_spe(), method = "spotlight",
    cell_type_col = "celltype", markers = markers, maxIter = 123)
  expect_identical(out, "model")
  expect_equal(received$mgs$gene, c("G1", "G2"))
  expect_identical(received$mgs$cluster, c("A", "B"))
  expect_equal(received$mgs$mean.AUC, c(.9, .95))
  expect_identical(received$weight_id, "mean.AUC")
  expect_equal(received$maxIter, 123)
})

test_that("real deconvolution: spotlight", {
  local_test_state()
  withr::local_dir(withr::local_tempdir())
  local_test_device()
  fixture <- spotlight_fixture()
  # First-seen cell types differ from alphabetical order. This produces a
  # non-symmetric topic reference and exposes incorrect NNLS orientation.
  fixture$sce <- fixture$sce[, c(21:30, 1:20)]
  model <- build_model(fixture$sce, spatial_obj = fixture$spe,
    cell_type_col = "celltype", method = "spotlight", markers = spotlight_markers(), maxIter = 200)
  expect_true(all(c("mod", "topic") %in% names(model)))
  result <- deconvolute(fixture$spe, signature = model, method = "spotlight", return_object = FALSE)
  expect_method_result(result, fixture$spe, "spotlight")
  fractions <- result[rownames(fixture$truth), paste0("spotlight_", colnames(fixture$truth)), drop = FALSE]
  expect_gte(mean(max.col(fractions) == max.col(fixture$truth)), .8)
  reordered <- model
  reordered$topic <- model$topic[c(3, 1, 2), c(2, 3, 1), drop = FALSE]
  again <- deconvolute(fixture$spe, signature = reordered, method = "spotlight",
                      return_object = FALSE)
  expect_equal(again[, colnames(result), drop = FALSE], result, tolerance = 1e-8)
})
