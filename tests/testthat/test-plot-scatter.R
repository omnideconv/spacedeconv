test_that("scatter labels can describe arbitrary numeric features or be omitted", {
  spe <- result_spe()
  spe$TF_activity <- c(-3, -1, 0, 2, 4, 3)
  p <- plot_scatter(spe, "TF_activity", "quantiseq_A", xlab = "TF activity",
                    ylab = "B cells", title = NULL)
  expect_identical(p$labels$x, "TF activity")
  expect_identical(p$labels$y, "B cells")
  expect_null(p$labels$title)
  expect_renderable(p)
  default <- plot_scatter(spe, "TF_activity", "quantiseq_A")
  expect_identical(default$labels$x, "TF_activity")
  expect_identical(default$labels$y, "quantiseq_A")
  blank <- plot_scatter(spe, "TF_activity", "quantiseq_A", xlab = NULL, ylab = NULL)
  expect_null(blank$labels$x)
  expect_null(blank$labels$y)
})

test_that("scatter coordinate limits preserve user ranges and negative activity values", {
  spe <- result_spe()
  spe$activity1 <- c(-3, -1, 0, 2, 4, 3)
  spe$activity2 <- c(-2, 0, 1, 3, 2, 5)
  p <- plot_scatter(spe, "activity1", "activity2", fix_coords = TRUE)
  expect_equal(p$coordinates$limits$x, c(-3, 5))
  expect_equal(p$coordinates$limits$y, c(-3, 5))
  expect_equal(p$coordinates$ratio, 1)
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$data[[1]]$x, spe$activity1)
  expect_equal(built$data[[1]]$y, spe$activity2)
  expect_renderable(p)
  for (fixed in c(FALSE, TRUE)) {
    p <- plot_scatter(spe, "activity1", "activity2", fix_coords = fixed,
                      coord_range = c(-10, 10))
    expect_equal(p$coordinates$limits$x, c(-10, 10))
    expect_equal(p$coordinates$limits$y, c(-10, 10))
  }
  expect_error(plot_scatter(spe, "activity1", "activity2", coord_range = c(1, -1)),
               "increasing")
  p <- plot_scatter(spe, "quantiseq_A", "quantiseq_B", log_scale = TRUE,
                    fix_coords = TRUE, coord_range = c(.01, 1))
  expect_equal(p$coordinates$limits$x, c(.01, 1))
  expect_renderable(p)
  expect_error(plot_scatter(spe, "quantiseq_A", "quantiseq_B", log_scale = TRUE,
                            coord_range = c(0, 1)), "positive")
})
