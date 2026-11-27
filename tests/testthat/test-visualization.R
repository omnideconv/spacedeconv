test_that("spatial plots use the expected values and can be rendered", {
  local_test_state()
  spe <- result_spe()
  cases <- list(
    list(plot_spatial, list(result = "quantiseq_A"), "quantiseq_A", spe$quantiseq_A),
    list(plot_celltype, list(cell_type = "quantiseq_A"), "quantiseq_A", spe$quantiseq_A),
    list(plot_gene, list(gene = "G1"), "gene", as.numeric(assay(spe)[1, ])),
    list(plot_umi_count, list(), "nUMI", unname(colSums(assay(spe)))),
    list(plot_ndetected_genes, list(), "ndetected_genes", unname(colSums(assay(spe) > 0))),
    list(plot_comparison, list(cell_type_1 = "quantiseq_A", cell_type_2 = "quantiseq_B"),
         "comparison", log((spe$quantiseq_A + 1) / (spe$quantiseq_B + 1))))
  for (case in cases) {
    p <- do.call(case[[1]], c(list(spe = spe, density = FALSE), case[[2]]))
    expect_renderable(p)
    expect_equal(as.numeric(p$layers[[1]]$data[[case[[3]]]]), case[[4]])
  }
  p <- plot_most_abundant(spe, method = "quantiseq")
  expect_renderable(p)
  expect_equal(as.character(p$layers[[1]]$data$mostAbundant), rep(c("quantiseq_A", "quantiseq_B"), each = 3))
  p <- plot_overview(spe)
  expect_s3_class(p, "plotly")
  built <- plotly::plotly_build(p)
  expect_equal(as.numeric(built$x$data[[1]]$x), unname(spatialCoords(spe)[, 1]))
  expect_equal(as.numeric(built$x$data[[1]]$marker$color), unname(colSums(assay(spe))))
  expect_error(plot_spatial(spe, result = "absent"), "No result columns")
})

test_that("scatter and signature comparisons align by IDs rather than position", {
  local_test_state()
  spe <- result_spe()
  p <- plot_scatter(spe1 = spe, spe2 = spe[, 6:1], value1 = "quantiseq_A", value2 = "quantiseq_B")
  expect_renderable(p)
  expect_equal(p$data$value1 + p$data$value2, rep(1, 6))
  sig <- matrix(1:8, 4, dimnames = list(paste0("G", 1:4), c("A", "B")))
  p <- compare_signatures(sig, sig[4:1, 2:1])
  expect_renderable(p)
  expect_equal(p$data$signature1, p$data$signature2)
})

test_that("image overlays, smoothing, transforms and method panels render", {
  local_test_state()
  spe <- result_spe()
  expect_renderable(plot_celltype(spe, "quantiseq_A", show_image = TRUE,
    smooth = TRUE, density = FALSE))
  p <- plot_gene(spe, "G1", transform_scale = "log", density = FALSE)
  expect_equal(p$layers[[1]]$data$gene, log1p(as.numeric(assay(spe)[1, ])))
  expect_renderable(p)
  expect_renderable(plot_spatial(spe, "quantiseq", density = FALSE))
})
