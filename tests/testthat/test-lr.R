spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))

spe <- preprocess(spe, remove_mito = TRUE)


test_that("LR works", {
  res <- get_lr(spe)

  expect_true(any(grepl("lr_", colnames(colData(res)))))
})
