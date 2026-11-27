test_that("configured Python environment provides required modules", {
  config <- reticulate::py_config()
  expect_true(reticulate::py_available())
  envname <- getOption("spacedeconv.conda_env", "spacedeconv-env")
  expect_true(grepl(envname, config$python, fixed = TRUE))
  expect_true(reticulate::py_module_available("igraph"))
  expect_true(reticulate::py_module_available("leidenalg"))
  expect_true(reticulate::py_module_available("community")) # python louvain
  expect_true(reticulate::py_module_available("sklearn"))
  expect_true(reticulate::py_module_available("scanpy"))
  expect_true(reticulate::py_module_available("rectanglepy"))
})
