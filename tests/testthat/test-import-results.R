external_results <- function() {
  data.frame(T_cells = c(.9, .8, .7, .1, .2, .3),
    Macrophages = c(.1, .2, .3, .9, .8, .7), row.names = paste0("spot", 1:6))
}

test_that("external results align by spot ID and preserve the spatial object", {
  spe <- result_spe()
  results <- external_results() * 10 # Counts must not be normalized to fractions.
  for (input in list(results, as.matrix(results))) {
    out <- import_deconvolution_results(spe, input[c(4, 1, 6, 2, 5, 3), ], "omnideconv_bayesprism")
    expect_preserved(out, spe)
    expect_identical(SpatialExperiment::imgData(out), SpatialExperiment::imgData(spe))
    expect_identical(out$omnideconv_bayesprism_T_cells, results$T_cells)
    expect_identical(out$omnideconv_bayesprism_Macrophages, results$Macrophages)
    expect_setequal(available_results(out, "omnideconv_bayesprism"),
      c("omnideconv_bayesprism_T_cells", "omnideconv_bayesprism_Macrophages"))
  }
  expect_identical(names(colData(spe)), names(colData(result_spe())))
  single <- import_deconvolution_results(spe[, 1], results[1, 1, drop = FALSE], "external")
  expect_equal(single$external_T_cells, 9)
})

test_that("external imports reject invalid spot IDs and incomplete results", {
  spe <- small_spe()
  results <- as.matrix(external_results())
  for (ids in list(NULL, rep("spot1", 6), c("unknown", paste0("spot", 2:6)),
                   c(NA, paste0("spot", 2:6)), c("", paste0("spot", 2:6)))) {
    input <- results
    rownames(input) <- ids
    expect_error(import_deconvolution_results(spe, input, "external"), "spot IDs")
  }
  expect_error(import_deconvolution_results(spe, results[-1, ], "external"), "spot IDs")
  expect_error(import_deconvolution_results(spe, t(results), "external"), "spot IDs")
  colnames(spe)[2] <- colnames(spe)[1]
  expect_error(import_deconvolution_results(spe, results, "external"), "spot IDs")
})

test_that("external imports validate values, names and column collisions", {
  spe <- small_spe()
  results <- as.matrix(external_results())
  expect_error(import_deconvolution_results(NULL, results, "external"), "SpatialExperiment")
  for (input in list(NULL, list(A = 1:6), results[FALSE, ], results[, FALSE], 1:6)) {
    expect_error(import_deconvolution_results(spe, input, "external"), "non-empty matrix")
  }
  for (value in list(NA_real_, NaN, Inf, -Inf, -.1, "text")) {
    input <- results
    input[1, 1] <- value
    expect_error(import_deconvolution_results(spe, input, "external"), "finite.*numeric")
  }
  for (prefix in list(NULL, NA_character_, "", "two words", "1method", "expression", "cluster", c("a", "b"))) {
    expect_error(import_deconvolution_results(spe, results, prefix), "result_name")
  }
  for (celltypes in list(NULL, c("A", "A"), c("A B", "A.B"), c("A", ""), c("A", NA))) {
    input <- results
    colnames(input) <- celltypes
    expect_error(import_deconvolution_results(spe, input, "external"), "Cell-type")
  }
  colnames(results)[1] <- "T cells"
  out <- import_deconvolution_results(spe, results, "external")
  expect_equal(out$external_T.cells, results[, 1], ignore_attr = TRUE)
  expect_error(import_deconvolution_results(out, results, "external"), "already exist")
  expect_equal(out$external_T.cells, results[, 1], ignore_attr = TRUE)
  second <- import_deconvolution_results(out, results, "another")
  expect_preserved(second, out)
  expect_equal(second$another_T.cells, out$external_T.cells)
})

test_that("custom prefix selection respects name boundaries", {
  spe <- import_deconvolution_results(small_spe(), external_results(), "external")
  spe <- import_deconvolution_results(spe, external_results(), "external2")
  expect_identical(available_results(spe, "external"), c("external_Macrophages", "external_T_cells"))
  expect_identical(available_results(spe, "external_"), available_results(spe, "external"))
  expect_length(available_results(spe, "extern"), 0)
  expect_error(cluster(spe, spmethod = "extern"), "No result columns")
  spe$invalid_score <- letters[1:6]
  expect_error(cluster(spe, spmethod = "invalid"), "finite numeric")
})

test_that("imported results support clustering and inferred feature prefixes", {
  local_test_state()
  spe <- import_deconvolution_results(small_spe(), external_results()[6:1, ], "omnideconv_bayesprism")
  for (method in c("kmeans", "hclust")) {
    out <- cluster(spe, method = method, spmethod = "omnideconv_bayesprism",
      nclusters = 2, dist_method = "euclidean")
    clusterid <- "cluster_omnideconv_bayesprism_nclusters_2"
    labels <- out[[clusterid]]
    expect_length(unique(labels[1:3]), 1)
    expect_length(unique(labels[4:6]), 1)
    expect_false(labels[1] == labels[4])
    expect_preserved(out, spe)
    features <- get_cluster_features(out, clusterid, topn = 1)
    expect_identical(names(features[[as.character(labels[1])]]), "omnideconv_bayesprism_T_cells")
    expect_identical(names(features[[as.character(labels[4])]]), "omnideconv_bayesprism_Macrophages")
  }
})

test_that("imported results support individual and grouped spatial plots", {
  local_test_state()
  spe <- import_deconvolution_results(small_spe(), external_results()[6:1, ], "omnideconv_bayesprism")
  columns <- available_results(spe, "omnideconv_bayesprism")
  p <- plot_spatial(spe, columns[1], density = FALSE)
  expect_renderable(p)
  expect_equal(p$layers[[1]]$data[[columns[1]]], spe[[columns[1]]])
  for (selection in list("omnideconv_bayesprism", columns)) {
    panels <- plot_spatial(spe, selection, density = FALSE)
    expect_s3_class(panels, "patchwork")
    expect_no_error(patchwork::patchworkGrob(panels))
    for (i in seq_along(columns)) {
      expect_equal(panels[[i]]$layers[[1]]$data[[columns[i]]], spe[[columns[i]]])
    }
  }
  p <- plot_most_abundant(spe, method = "omnideconv_bayesprism")
  expect_renderable(p)
  expect_equal(as.character(p$layers[[1]]$data$mostAbundant),
    rep(c("omnideconv_bayesprism_T_cells", "omnideconv_bayesprism_Macrophages"), each = 3))
  one <- import_deconvolution_results(small_spe(), external_results()[, 1, drop = FALSE], "single")
  expect_renderable(plot_spatial(one, "single", density = FALSE))
})
