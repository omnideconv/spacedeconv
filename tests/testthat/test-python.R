library(spacedeconv)

test_that("Python environment creation works", {
  expect_true(reticulate::py_available())
  expect_true(grepl("r-spacedeconv", reticulate::py_config()$python))
  expect_true(reticulate::py_module_available("scaden"))
  expect_true(reticulate::py_module_available("igraph"))
  expect_true(reticulate::py_module_available("leidenalg"))
  expect_true(reticulate::py_module_available("community")) # python louvain
  expect_true(reticulate::py_module_available("sklearn"))
  #expect_true(reticulate::py_module_available("python.app"))
  expect_true(reticulate::py_module_available("scanpy"))
  expect_true(reticulate::py_module_available("cell2location"))
  expect_true(reticulate::py_module_available("metacells"))
  expect_true(reticulate::py_module_available("autogenes"))
})
