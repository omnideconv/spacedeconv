test_that("DecoupleR Reference can be computed (progeny)", {
  ref <- get_decoupleR_reference(method = "progeny", organism = "human")
  expect_contains(colnames(ref), "source")
  expect_contains(colnames(ref), "target")
})


spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))

test_that("DecoupleR quantification works (progeny)", {
  ref <- get_decoupleR_reference(method = "progeny", organism = "human")
  spe <- preprocess(spe, remove_mito = TRUE)
  spe <- normalize(spe)
  res <- compute_activities(spe, reference = ref)

  expect_s4_class(res, "SpatialExperiment")
  expect_contains(colnames(colData(res)), "progeny_WNT")
})
