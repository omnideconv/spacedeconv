test_that("Python environment creation works", {
  init_python()
  expect_identical(info="Python available", object = reticulate::py_available(), expected = TRUE)
})
