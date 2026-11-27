test_that("annotations and available results preserve data and order", {
  spe <- result_spe()
  expect_equal(available_results(spe, "quantiseq"), c("quantiseq_A", "quantiseq_B"))
  expect_length(available_results(spe, "absent"), 0)
  out <- annotate_spots(spe, c("spot5", "spot2"), name = "selected")
  expect_identical(out$selected, c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE))
  expect_preserved(out, spe)
  out <- addCustomAnnotation(spe, "note", letters[1:6])
  expect_identical(out$note, letters[1:6])
  expect_preserved(out, spe)
  expect_error(addCustomAnnotation(spe, "note", 1:2), "length")
})

test_that("aggregation and cell count scaling have known numerical answers", {
  spe <- result_spe()
  out <- aggregate_results(spe, c("quantiseq_A", "quantiseq_B"), name = "total")
  expect_equal(out$total, rep(1, 6))
  expect_preserved(out, spe)
  removed <- aggregate_results(spe, c("quantiseq_A", "quantiseq_B"), remove = TRUE)
  expect_false(any(c("quantiseq_A", "quantiseq_B") %in% names(colData(removed))))
  expect_warning(legacy <- aggregate_results(spe, cell_type_1 = "quantiseq_A",
    cell_type_2 = "quantiseq_B"), "deprecated")
  expect_equal(legacy$quantiseq_A_quantiseq_B, rep(1, 6))
  expect_error(aggregate_results(spe, c("quantiseq_A", "quantiseq_A")), "Duplicate")
  expect_error(aggregate_results(spe, c("quantiseq_A", "missing")), "not available")
  spe$n_cells <- 1:6
  out <- scale_cell_counts(spe, "quantiseq_A", "n_cells", resName = "absolute")
  expect_equal(out$absolute, spe$quantiseq_A * (1:6))
  expect_preserved(out, spe)
  expect_error(scale_cell_counts(spe, "missing", "n_cells"), "does not exist")
  expect_error(scale_cell_counts(spe, "quantiseq_A", "missing"), "not available")
})

test_that("all normalization formulas preserve counts and give known values", {
  spe <- small_spe()
  x <- as.matrix(assay(spe))
  cpm <- sweep(x, 2, colSums(x), "/") * 1e6
  expected <- list(cpm = cpm, logcpm = log1p(cpm),
    logpf = log(sweep(x + 1, 1, rowMeans(x + 1) + 1, "/")))
  for (method in names(expected)) {
    out <- normalize(spe, method = method)
    expect_equal(as.matrix(assay(out, method)), expected[[method]], tolerance = 1e-10)
    expect_preserved(out, spe)
  }
  expect_error(normalize(NULL), "null")
  expect_error(normalize(spe, method = "unknown"), "arg")
  expect_error(normalize(spe, assay = "absent"))
})

test_that("spatial subsetting handles boundaries, one spot and empty selections", {
  spe <- small_spe()
  expect_equal(subsetSPE(spe), spe)
  out <- subsetSPE(spe, colRange = c(10, 0), rowRange = c(0, 10))
  expect_identical(colnames(out), c("spot1", "spot2", "spot4", "spot5"))
  expect_equal(spatialCoords(out), spatialCoords(spe)[c(1, 2, 4, 5), , drop = FALSE])
  expect_identical(colnames(subsetSPE(spe, c(0, 0), c(0, 0))), "spot1")
  expect_equal(ncol(subsetSPE(spe, c(100, 200))), 0)
  spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = assay(spe)),
    spatialCoords = spatialCoords(spe), sample_id = rep(c("one", "two"), each = 3))
  expect_equal(filter_sample_id(spe, "two"), spe[, 4:6])
  expect_error(filter_sample_id(spe, NULL))
  expect_equal(filter_sample_id(small_spe(), NULL), small_spe())
})

test_that("single cell subsetting is reproducible and respects both strategies", {
  local_test_state()
  sce <- small_sce()
  for (scenario in c("even", "mirror")) {
    out <- subsetSCE(sce, "celltype", scenario, ncells = 4, seed = 7)
    expect_equal(as.integer(table(out$celltype)), c(2L, 2L))
    expect_equal(out, subsetSCE(sce, "celltype", scenario, ncells = 4, seed = 7))
    expect_equal(assay(out), assay(sce)[, colnames(out), drop = FALSE])
  }
  expect_equal(subsetSCE(sce, "celltype", ncells = 20, notEnough = "asis"), sce)
  expect_equal(ncol(subsetSCE(sce, "celltype", ncells = 20, notEnough = "remove")), 0)
  expect_error(subsetSCE(sce, "missing"), "can't be found")
  expect_error(subsetSCE(sce, "celltype", seed = "bad"), "numeric")
})

test_that("dataset summaries validate inputs and describe real dimensions", {
  spe <- small_spe()
  expect_no_error(print_info(sce = small_sce(), spe = spe,
    signature = as.matrix(assay(spe))[, 1:2]))
  expect_error(print_info(sce = 1), "datatype")
  expect_error(print_info(spe = 1), "SpatialExperiment")
})
