library(spacedeconv)

test_that("Python environment creation works", {
  expect_true(reticulate::py_available())
  envname <- getOption("omnideconv.conda_env", "spacedeconv")
  expect_true(grepl(envname, reticulate::py_config()$python, fixed = TRUE))
  expect_true(reticulate::py_module_available("igraph"))
  expect_true(reticulate::py_module_available("leidenalg"))
  expect_true(reticulate::py_module_available("community")) # python louvain
  expect_true(reticulate::py_module_available("sklearn"))
  expect_true(reticulate::py_module_available("scanpy"))
})
