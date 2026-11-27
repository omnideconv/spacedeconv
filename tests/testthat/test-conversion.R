test_that("human and mouse ortholog conversion maps genes and complexes", {
  human <- c("BRCA1", "TP53", "BRCA1_TP53")
  mouse <- c("Brca1", "Trp53", "Brca1_Trp53")
  expect_equal(convert_human_to_mouse(human)$Mouse_symbol, mouse)
  expect_equal(convert_mouse_to_human(mouse)$Human_symbol, human)
})

test_that("AnnData conversion preserves values, names and metadata", {
  spe <- small_spe()
  ad <- spe_to_ad(spe)
  expect_s3_class(ad, "AnnDataR6")
  expect_equal(as.matrix(ad$X), t(as.matrix(assay(spe))), ignore_attr = TRUE)
  expect_equal(as.character(ad$obs_names), colnames(spe))
  expect_equal(as.character(ad$var_names), rownames(spe))
  sce <- anndata_to_singlecellexperiment(ad)
  expect_equal(as.matrix(assay(sce, "X")), as.matrix(assay(spe)))
  expect_equal(as.character(sce$sample_id), as.character(spe$sample_id))
  # Supply the standard spatial AnnData fields (spe_to_ad only exports expression).
  obs <- as.data.frame(colData(spe))
  obs$sample <- factor(rep("sample01", ncol(spe)))
  ad <- anndata::AnnData(X = Matrix::t(assay(spe)), obs = obs,
    var = data.frame(row.names = rownames(spe)),
    obsm = list(spatial = spatialCoords(spe)),
    uns = list(spatial = list(sample01 = list(images = list(lowres = array(1, c(30, 30, 3))),
      scalefactors = list(tissue_lowres_scalef = 1)))))
  out <- anndata_to_spatialexperiment(ad)
  expect_equal(as.matrix(assay(out)), as.matrix(assay(spe)))
  expect_equal(unname(spatialCoords(out)), unname(spatialCoords(spe)))
  expect_equal(SpatialExperiment::scaleFactors(out), 1)
  expect_error(spe_to_ad(NULL), "missing")
  expect_error(spe_to_ad(spe, "absent"), "not available")
})

test_that("Seurat spatial conversion preserves counts, coordinates and image", {
  local_test_state()
  spe <- small_spe()
  obj <- Seurat::CreateSeuratObject(counts = assay(spe), assay = "Spatial")
  coords <- data.frame(tissue = 1, row = 1:6, col = 1:6,
    imagerow = spatialCoords(spe)[, 2], imagecol = spatialCoords(spe)[, 1],
    row.names = colnames(spe))
  obj[["slice1"]] <- methods::new("VisiumV1", image = array(1, c(30, 30, 3)),
    scale.factors = Seurat::scalefactors(spot = 1, fiducial = 1, hires = 1, lowres = 1),
    coordinates = coords, spot.radius = 1, assay = "Spatial", key = "slice1_")
  out <- seurat_to_spatialexperiment(obj)
  expect_equal(assay(out), assay(spe))
  expect_equal(spatialCoords(out)[, "pxl_col_in_fullres"], spatialCoords(spe)[, "pxl_col_in_fullres"])
  expect_equal(spatialCoords(out)[, "pxl_row_in_fullres"], spatialCoords(spe)[, "pxl_row_in_fullres"])
  expect_equal(SpatialExperiment::scaleFactors(out), 1)
})
