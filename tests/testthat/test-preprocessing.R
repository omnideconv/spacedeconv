test_that("preprocessing selects exact spots and genes without changing input", {
  spe <- small_spe()
  rownames(spe) <- c("MT-A", "G2", "G3", "G4")
  original <- spe
  totals <- colSums(assay(spe))
  out <- preprocess(spe, min_umi = 8, max_umi = 12, remove_mito = TRUE)
  expect_identical(colnames(out), colnames(spe)[totals >= 8 & totals <= 12])
  expect_equal(assay(out), assay(spe)[2:4, totals >= 8 & totals <= 12, drop = FALSE])
  expect_equal(spatialCoords(out), spatialCoords(spe)[totals >= 8 & totals <= 12, , drop = FALSE])
  expect_identical(spe, original)
  expect_warning(preprocess(spe, min_umi = 0), "mitochondrial")
  expect_error(preprocess(NULL), "provide an object")
})

test_that("preprocessing removes zero genes and keeps the strongest duplicate", {
  spe <- small_spe()
  assay(spe)[1, ] <- 0
  rownames(spe) <- c("zero", "duplicate", "duplicate", "G4")
  out <- preprocess(spe, min_umi = 0)
  expect_setequal(rownames(out), c("duplicate", "G4"))
  expect_equal(as.numeric(assay(out)["duplicate", ]), as.numeric(assay(spe)[2, ]))
  expect_s4_class(preprocess(small_sce(), min_umi = 0), "SingleCellExperiment")
})
