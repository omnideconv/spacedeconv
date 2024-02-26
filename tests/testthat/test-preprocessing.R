sce <- readRDS(system.file("testdata", "sce.rds", package = "spacedeconv"))
spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))

test_that("Function rejects invalid inputs", {
  expect_error(spacedeconv::preprocess(NULL), "Please provide an object")
})


test_that("UMI threshold filtering works", {
  filtered_sce <- preprocess(sce, min_umi = 500, max_umi = 1500, remove_mito = TRUE)
  filtered_spe <- preprocess(spe, min_umi = 500, max_umi = 1500, remove_mito = TRUE)

  expect_true(all(colSums(assay(filtered_sce)) >= 500))
  expect_true(all(colSums(assay(filtered_sce)) <= 1500))
  expect_true(all(colSums(assay(filtered_spe)) >= 500))
  expect_true(all(colSums(assay(filtered_spe)) <= 1500))
})


test_that("Mitochondrial genes are removed when requested", {
  processed_sce <- preprocess(sce, remove_mito = TRUE)
  processed_spe <- preprocess(spe, remove_mito = TRUE)

  expect_false(any(grepl("^MT-", rownames(assay(processed_sce)))))
  expect_false(any(grepl("^MT-", rownames(assay(processed_spe)))))
})
