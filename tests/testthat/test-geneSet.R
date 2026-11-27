test_that("gene set scores equal the sum of log1p expression", {
  spe <- small_spe()
  out <- gene_set_score(spe, c("G1", "G3", "missing"), assay = "counts", name = "score")
  expect_equal(unname(out$score), unname(colSums(log1p(as.matrix(assay(spe))[c(1, 3), ]))))
  expect_preserved(out, spe)
  expect_error(gene_set_score(spe, "absent"), "No gene")
  expect_error(gene_set_score(spe, NULL), "non-null")
  expect_error(gene_set_score(NULL, "G1"), "SpatialExperiment")
})

test_that("ligand receptor min and product scores have known values", {
  spe <- normalize(small_spe())
  resource <- data.frame(source_genesymbol = c("G1", "G2"), target_genesymbol = c("G3", "G4"))
  x <- as.matrix(assay(spe, "cpm"))
  for (method in c("min", "product")) {
    out <- get_lr(spe, resource = resource, method = method)
    expected <- if (method == "min") pmin(x[1, ], x[3, ]) else x[1, ] * x[3, ]
    expect_equal(out$lr_G1.G3, unname(expected))
    expect_preserved(out, spe)
  }
  expect_error(get_lr(small_spe(), resource), "cpm")
  expect_error(get_lr(spe, resource, method = "unknown"), "not supported")
})
