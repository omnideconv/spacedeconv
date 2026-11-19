test_that("spatial correlations preserve explicit order and support arbitrary numeric columns", {
  local_test_device()
  spe <- result_spe()
  spe$activity <- c(-3, 1, 0, 8, 2, -1)
  columns <- c("activity", "quantiseq_B", "quantiseq_A")
  categories <- c(quantiseq_A = "Cells", activity = "Activity", quantiseq_B = "Cells")
  labels <- c(quantiseq_B = "B", quantiseq_A = "A", activity = "TF")
  for (method in c("pearson", "spearman", "kendall")) {
    out <- spatialcorr(spe, columns, cor_method = method, categories = categories,
      labels = labels, category_colors = c(Cells = "blue", Activity = "orange"),
      row_names_gp = grid::gpar(fontsize = 8))
    expect_equal(out$corr, cor(as.matrix(colData(spe)[, columns]), method = method))
    expect_equal(ComplexHeatmap::row_order(out$heatmap), 1:3)
    expect_equal(ComplexHeatmap::column_order(out$heatmap), 1:3)
    expect_identical(out$heatmap@ht_list[[1]]@row_names_param$labels, c("TF", "B", "A"))
    expect_identical(out$used_spots, colnames(spe))
    expect_length(out$excluded_spots, 0)
    expect_length(out$excluded_columns, 0)
  }
  out <- spatialcorr(spe, columns)
  expect_s4_class(out$heatmap, "HeatmapList")
  expect_null(out$padj)
})

test_that("spatial correlations report exclusions and use the same spots for every pair", {
  local_test_device()
  spe <- result_spe()
  spe$activity <- c(NA, Inf, 2, 4, 1, 3)
  spe$constant <- rep(0, 6)
  columns <- c("quantiseq_A", "activity", "quantiseq_B", "constant")
  warnings <- character()
  out <- withCallingHandlers(spatialcorr(spe, columns), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(warnings, 2)
  expect_match(warnings[1], "2 spots")
  expect_match(warnings[2], "constant")
  expect_identical(out$excluded_spots, colnames(spe)[1:2])
  expect_identical(out$used_spots, colnames(spe)[3:6])
  expect_identical(out$excluded_columns, "constant")
  expect_equal(out$corr, cor(as.matrix(colData(spe)[3:6, columns[1:3]])))
  spe$activity <- c(1, rep(NA_real_, 5))
  expect_error(suppressWarnings(spatialcorr(spe, columns)), "two complete spots")
  expect_error(suppressWarnings(spatialcorr(spe, c("quantiseq_A", "constant"))), "two nonconstant")
})

test_that("spatial correlations reject invalid selections and ambiguous annotations", {
  spe <- result_spe()
  columns <- c("quantiseq_A", "quantiseq_B")
  expect_error(spatialcorr(spe, c("missing", "quantiseq_A")), "do not exist")
  expect_error(spatialcorr(spe, rep("quantiseq_A", 2)), "distinct")
  expect_error(spatialcorr(spe, c("sample_id", "quantiseq_A")), "numeric")
  expect_error(spatialcorr(spe, columns, categories = c("A", "B")), "named")
  expect_error(spatialcorr(spe, columns, labels = c(quantiseq_A = "A")), "covering")
  expect_error(spatialcorr(spe, columns, cluster_rows = TRUE), "Cannot override")
  expect_error(spatialcorr(spe, columns, cor_method = "bad"), "arg")
  expect_error(spatialcorr(spe, columns, category_colors = c(A = "red")), "requires categories")
})

test_that("expression correlation passes the correct values to its plot", {
  local_test_device()
  sig <- cbind(A = c(1, 4, 2, 8, 5, 3), B = c(6, 2, 4, 1, 3, 5), C = c(2, 3, 6, 4, 1, 8))
  captured <- NULL
  local_mocked_bindings(corrplot = function(corr, ...) { captured <<- corr })
  for (method in c("pearson", "spearman")) {
    corr_expr(sig, cor_method = method)
    expect_equal(captured, cor(sig, method = method))
  }
  corr_expr(sig, log = TRUE)
  expect_equal(captured, cor(log1p(sig)))
  expect_error(corr_expr(sig, cor_method = "bad"), "cor_method")
})

test_that("custom row-label headings render and survive redrawing", {
  local_test_device()
  spe <- result_spe()
  columns <- c("quantiseq_A", "quantiseq_B")
  out <- spatialcorr(spe, columns, labels = c(quantiseq_A = "First", quantiseq_B = "Second"),
    categories = c(quantiseq_A = "A", quantiseq_B = "B"), label_title = "Feature | Description")
  expect_true(out$heatmap@ht_list[[1]]@right_annotation@anno_list[[1]]@name_param$show)
  expect_no_error(ComplexHeatmap::draw(out$heatmap))
  expect_equal(ComplexHeatmap::row_order(out$heatmap), 1:2)
  expect_equal(out$corr, cor(as.matrix(colData(spe)[, columns])))
  expect_error(spatialcorr(spe, columns, label_title = c("A", "B")), "single string")
})

test_that("row labels and their heading can be hidden independently", {
  local_test_device()
  spe <- result_spe()
  columns <- c("quantiseq_A", "quantiseq_B")
  no_heading <- spatialcorr(spe, columns, label_title = NULL)
  expect_null(no_heading$heatmap@ht_list[[1]]@right_annotation)
  expect_true(no_heading$heatmap@ht_list[[1]]@row_names_param$show)
  hidden <- spatialcorr(spe, columns, show_row_names = FALSE)
  expect_null(hidden$heatmap@ht_list[[1]]@right_annotation)
  expect_false(hidden$heatmap@ht_list[[1]]@row_names_param$show)
  expect_equal(hidden$corr, no_heading$corr)
})

test_that("custom legend and heading styles survive moving legends below the heatmap", {
  local_test_device()
  spe <- result_spe()
  columns <- c("quantiseq_A", "quantiseq_B")
  out <- spatialcorr(spe, columns,
    categories = c(quantiseq_A = "A", quantiseq_B = "B"),
    annotation_legend_param = list(nrow = 2, labels_gp = grid::gpar(fontsize = 20)),
    label_title_gp = grid::gpar(fontsize = 20, fontface = "bold"))
  redrawn <- ComplexHeatmap::draw(out$heatmap@ht_list[[1]],
    heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  expect_identical(redrawn@heatmap_legend_param$side, "bottom")
  expect_identical(redrawn@annotation_legend_param$side, "bottom")
  ht <- redrawn@ht_list[[1]]
  legend <- ht@top_annotation@anno_list$Category@legend_param
  expect_equal(legend$nrow, 2)
  expect_equal(legend$labels_gp$fontsize, 20)
  heading <- ht@right_annotation@anno_list$Labels@name_param
  expect_equal(heading$gp$fontsize, 20)
  expect_identical(heading$just, c("left", "bottom"))
  expect_equal(heading$x, grid::unit(0, "npc"))
  expect_equal(out$corr, cor(as.matrix(colData(spe)[, columns])))
  expect_error(spatialcorr(spe, columns, annotation_legend_param = list(2)), "named list")
  expect_error(spatialcorr(spe, columns, label_title_gp = 20), "gpar")
})


test_that("category frames survive redrawing and follow retained contiguous blocks", {
  local_test_device()
  spe <- result_spe()
  spe$activity <- c(-3, 1, 0, 8, 2, -1)
  spe$constant <- rep(0, 6)
  columns <- c("quantiseq_A", "constant", "quantiseq_B", "activity")
  categories <- setNames(c("Cells", "Removed", "Cells", "Activity"), columns)
  expect_warning(out <- spatialcorr(spe, columns, categories = categories,
    draw_param = list(heatmap_legend_side = "bottom", annotation_legend_side = "bottom")),
    "constant")
  expect_identical(out$heatmap@heatmap_legend_param$side, "bottom")
  expect_identical(out$heatmap@annotation_legend_param$side, "bottom")
  for (i in 1:2) {
    blocks <- grid::grid.get("category_blocks")
    expect_equal(as.numeric(blocks$width), c(2/3, 1/3))
    expect_equal(as.numeric(blocks$x), c(1/3, 5/6))
    expect_equal(blocks$gp$lwd, 3)
    ComplexHeatmap::draw(out$heatmap)
  }
  # A single category has a frame but no internal boundaries.
  expect_no_error(spatialcorr(spe, columns[c(1, 3)],
    categories = setNames(c("Cells", "Cells"), columns[c(1, 3)])))
  expect_error(spatialcorr(spe, columns, draw_param = list("bottom")), "named list")
  expect_error(spatialcorr(spe, columns, draw_param = list(cluster_rows = TRUE)), "cannot override")
})
