# Fresh, deterministic objects: no shared mutable S4/Python fixtures.
small_spe <- function() {
  x <- matrix(c(1, 3, 0, 2, 4, 2, 1, 0, 3, 1, 5, 2,
                2, 6, 2, 1, 8, 1, 3, 4, 3, 5, 1, 7), nrow = 4,
              dimnames = list(c("G1", "G2", "G3", "G4"), paste0("spot", 1:6)))
  coords <- cbind(pxl_col_in_fullres = c(0, 10, 20, 0, 10, 20),
                  pxl_row_in_fullres = c(0, 0, 0, 10, 10, 10))
  rownames(coords) <- colnames(x)
  image <- SpatialExperiment::SpatialImage(as.raster(matrix("white", 30, 30)))
  SpatialExperiment::SpatialExperiment(
    assays = list(counts = as(x, "dgCMatrix")), spatialCoords = coords,
    colData = S4Vectors::DataFrame(sample_id = rep("sample01", 6),
      in_tissue = rep(1L, 6), row.names = colnames(x)),
    imgData = S4Vectors::DataFrame(sample_id = "sample01", image_id = "lowres",
      data = I(list(image)), scaleFactor = 1))
}

small_sce <- function() {
  spe <- small_spe()
  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = SummarizedExperiment::assay(spe)),
    colData = S4Vectors::DataFrame(celltype = rep(c("A", "B"), each = 3),
      sample_id = rep("sample01", 6), row.names = colnames(spe)))
}

result_spe <- function() {
  spe <- small_spe()
  spe$quantiseq_A <- c(.9, .8, .7, .1, .2, .3)
  spe$quantiseq_B <- 1 - spe$quantiseq_A
  spe
}

expect_preserved <- function(actual, original) {
  expect_identical(SummarizedExperiment::assay(actual, "counts"),
                   SummarizedExperiment::assay(original, "counts"))
  expect_identical(colnames(actual), colnames(original))
  expect_identical(SpatialExperiment::spatialCoords(actual),
                   SpatialExperiment::spatialCoords(original))
  expect_equal(SummarizedExperiment::colData(actual)[, names(SummarizedExperiment::colData(original)), drop = FALSE],
               SummarizedExperiment::colData(original))
}

local_test_state <- function(seed = 42, .local_envir = parent.frame()) {
  withr::local_seed(seed, .local_envir = .local_envir)
  withr::local_options(list(lifecycle_verbosity = "quiet"), .local_envir = .local_envir)
}

local_test_device <- function(.local_envir = parent.frame()) {
  grDevices::pdf(file = NULL)
  device <- grDevices::dev.cur()
  withr::defer(grDevices::dev.off(device), envir = .local_envir)
}

expect_renderable <- function(p) {
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplotGrob(p))
}

local_python_state <- function(.local_envir = parent.frame()) {
  torch <- reticulate::import("torch", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  random <- reticulate::import("random", convert = FALSE)
  threads <- torch$get_num_threads()
  torch_state <- torch$get_rng_state()
  numpy_state <- np$random$get_state()
  random_state <- random$getstate()
  withr::defer({
    torch$set_num_threads(threads)
    torch$set_rng_state(torch_state)
    np$random$set_state(numpy_state)
    random$setstate(random_state)
  }, envir = .local_envir)
  torch$set_num_threads(1L)
  torch$manual_seed(42L)
  np$random$seed(42L)
  random$seed(42L)
}
