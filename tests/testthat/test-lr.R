spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))

# subset for speed
spe <- spe[1:10, ]

test_that("LR works", {
  spe <- get_lr(spe)

  expect_true(any(grepl("lr", colnames(colData(spe)))))
})
