test_that("result insertion aligns shuffled and incomplete results by spot ID", {
  spe <- small_spe()
  result <- matrix(1:12, 6, dimnames = list(rev(colnames(spe)), c("method_A", "method_B")))
  out <- addResultToObject(spe, result)
  expect_equal(unname(out$method_A), unname(result[colnames(spe), "method_A"]))
  expect_preserved(out, spe)
  partial <- result[c(1, 3, 6), , drop = FALSE]
  out <- addResultToObject(spe, partial)
  expect_equal(out$method_A, unname(partial[match(colnames(spe), rownames(partial)), "method_A"]), ignore_attr = TRUE)
  expect_identical(unname(which(is.na(out$method_A))), c(2L, 3L, 5L))
  expect_error(addResultToObject(spe, result[c(1, 1), , drop = FALSE]), "unique")
  rownames(result)[1] <- "unknown"
  expect_error(addResultToObject(spe, result), "unknown|Unknown")
})

test_that("immunedeconv conversion retains single spots and cell type names", {
  input <- data.frame(cell_type = c("A", "B"), spot2 = c(.2, .8), spot1 = c(.7, .3))
  out <- convertImmunedeconvMatrix(input)
  expect_equal(out, matrix(c(.2, .7, .8, .3), 2,
    dimnames = list(c("spot2", "spot1"), c("A", "B"))))
  expect_equal(dim(convertImmunedeconvMatrix(input[, 1:2])), c(1L, 2L))
  expect_equal(colnames(attachToken(out, "custom")), c("custom_A", "custom_B"))
  expect_equal(unname(attachToken(out, "custom")), unname(out))
})

test_that("internal data helpers preserve values and identify unusable spots", {
  spe <- small_spe()
  dense <- spe
  assay(dense) <- as.matrix(assay(dense))
  expect_equal(check_datatype(dense), spe)
  assay(spe)[, 2] <- 0
  out <- removeZeroExpression(spe)
  expect_identical(colnames(out), colnames(spe)[-2])
  expect_equal(spatialCoords(out), spatialCoords(spe)[-2, , drop = FALSE])
  expect_false(checkRowColumn(spe))
  rownames(spe) <- NULL
  expect_true(checkRowColumn(spe))
  sce <- small_sce()
  matrices <- getMatricesFromSCE(sce, "celltype")
  expect_equal(matrices$counts, assay(sce))
  expect_equal(matrices$cell_type_annotation, sce$celltype)
  expect_true(checkENSEMBL(c("ENSG0000001", "G1")))
  expect_false(checkENSEMBL(c("G1", "G2")))
})

test_that("ligand receptor complexes handle zero expression per spot", {
  x <- data.frame(spot1 = c(0, 2, 4), spot2 = c(6, 2, 4),
    row.names = c("A", "B", "C"))
  expect_equal(as.numeric(get_ligand_expression("A_B~C", x)), c(0, 4))
  expect_equal(as.numeric(get_receptor_expression("C~A_B", x)), c(0, 4))
  expect_equal(as.numeric(get_ligand_expression("A_missing~C", x)), c(0, 0))
  expect_equal(as.numeric(get_receptor_expression("C~missing_B", x)), c(0, 0))
})
