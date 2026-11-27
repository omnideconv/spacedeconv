# All registered methods have a real test in the same default suite.
method_cases <- list(
  mcp_counter = list(), epic = list(), quantiseq = list(), xcell = list(),
  cibersort = list(licensed = TRUE), cibersort_abs = list(licensed = TRUE),
  timer = list(args = list(indications = rep("brca", 6))),
  consensus_tme = list(args = list(indications = rep("BRCA", 6))),
  abis = list(), estimate = list(),
  mmcp_counter = list(mouse = TRUE), seqimmucc = list(mouse = TRUE, args = list(algorithm = "LLSR")),
  dcq = list(mouse = TRUE), base = list(mouse = TRUE),
  rctd = list(reference = TRUE, args = list(n_cores = 1)),
  spotlight = list(reference = TRUE), card = list(reference = TRUE, package = "MuSiC"),
  spatialdwls = list(reference = TRUE),
  rectangle = list(reference = TRUE),
  cell2location = list(reference = TRUE), dot = list(reference = TRUE, package = "DOTr"))

reference_mixture <- function() {
  withr::local_seed(1701)
  ng <- 1200L
  nc <- 120L
  types <- rep(c("A", "B", "C"), each = nc / 3)
  profiles <- matrix(2, ng, 3, dimnames = list(paste0("Gene", seq_len(ng)), c("A", "B", "C")))
  for (i in 1:3) profiles[seq((i - 1) * 300 + 1, i * 300), i] <- 35
  profiles <- profiles * runif(ng, .5, 1.5)
  counts <- matrix(rpois(ng * nc, as.vector(profiles[, types])), ng,
    dimnames = list(rownames(profiles), paste0("cell", seq_len(nc))))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as(counts, "dgCMatrix")),
    colData = S4Vectors::DataFrame(celltype = types, sample_id = rep("sample01", nc), row.names = colnames(counts)))
  truth <- matrix(.1, 36, 3, dimnames = list(paste0("spot", 1:36), c("A", "B", "C")))
  truth[cbind(1:36, rep(1:3, each = 12))] <- .8
  means <- profiles %*% t(truth)
  sp <- matrix(rpois(length(means), as.vector(means)), ng, dimnames = dimnames(means))
  coords <- as.matrix(expand.grid(x = 1:6, y = 1:6))
  rownames(coords) <- colnames(sp)
  spe <- SpatialExperiment::SpatialExperiment(assays = list(counts = as(sp, "dgCMatrix")), spatialCoords = coords,
    sample_id = rep("sample01", ncol(sp)))
  list(sce = sce, spe = spe, truth = truth)
}

human_method_input <- function() {
  sce <- readRDS(system.file("testdata", "sce.rds", package = "spacedeconv", mustWork = TRUE))
  template <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv", mustWork = TRUE))[, 1:6]
  # Pool the reference instead of using sparse individual spots: TIMER needs
  # broad gene coverage and quanTIseq can be undefined for marker-poor spots.
  weights <- matrix(seq_len(ncol(sce) * 6) %% 11 + 1, ncol(sce), 6)
  bulk <- as.matrix(assay(sce) %*% weights)
  colnames(bulk) <- colnames(template)
  SpatialExperiment::SpatialExperiment(assays = list(counts = as(bulk, "dgCMatrix")),
    colData = colData(template), spatialCoords = spatialCoords(template),
    imgData = SpatialExperiment::imgData(template))
}

mouse_method_input <- function() {
  path <- system.file("extdata", "mouse_deconvolution", "sig_matr_seqImmuCC.txt",
    package = "immunedeconv", mustWork = TRUE)
  sig <- as.matrix(read.delim(path, row.names = 1, check.names = FALSE))
  # Real mouse symbols and profiles supplied by immunedeconv; no online orthology lookup.
  weights <- matrix(seq_len(ncol(sig) * 6) %% 7 + 1, ncol(sig), 6)
  weights <- sweep(weights, 2, colSums(weights), "/")
  x <- sig %*% weights
  colnames(x) <- paste0("mouse", 1:6)
  SpatialExperiment::SpatialExperiment(assays = list(counts = as(x, "dgCMatrix")),
    spatialCoords = cbind(x = 1:6, y = c(0, 1, 0, 1, 0, 1)), sample_id = rep("sample01", 6))
}

local_cibersort <- function(.local_envir = parent.frame()) {
  script <- Sys.getenv("SPACEDECONV_CIBERSORT_SCRIPT")
  signature <- Sys.getenv("SPACEDECONV_CIBERSORT_MATRIX")
  if (!file.exists(script) || !file.exists(signature)) {
    skip("CIBERSORT requires licensed files: SPACEDECONV_CIBERSORT_SCRIPT and SPACEDECONV_CIBERSORT_MATRIX")
  }
  local_config(.local_envir)
  set_cibersort_binary(script)
  set_cibersort_mat(signature)
}

local_config <- function(.local_envir = parent.frame()) {
  for (env in list(config_env, get("config_env", asNamespace("immunedeconv")))) {
    local({
      target <- env
      old <- as.list(target)
      withr::defer({rm(list = ls(target, all.names = TRUE), envir = target); list2env(old, target)},
        envir = .local_envir)
    })
  }
}

expect_method_result <- function(result, spe, method) {
  result <- as.matrix(result)
  expect_true(is.numeric(result))
  expect_equal(nrow(result), ncol(spe))
  expect_gt(ncol(result), 0)
  expect_setequal(rownames(result), colnames(spe))
  expect_false(anyDuplicated(rownames(result)) > 0)
  expect_false(anyDuplicated(colnames(result)) > 0)
  expect_true(all(startsWith(colnames(result), paste0(method, "_"))))
  expect_true(all(is.finite(result)))
  if (method %in% c("quantiseq", "epic", "cibersort", "rctd", "card", "spotlight", "spatialdwls", "cell2location", "dot", "seqimmucc")) {
    expect_true(all(result >= -1e-8))
  }
  if (method %in% c("quantiseq", "epic", "cibersort", "rctd", "card", "cell2location")) {
    expect_equal(unname(rowSums(result)), rep(1, nrow(result)), tolerance = 1e-4)
  }
  out <- addResultToObject(spe, result)
  expect_preserved(out, spe)
  expect_equal(as.matrix(colData(out)[, make.names(colnames(result), unique = TRUE), drop = FALSE]),
    result[colnames(spe), , drop = FALSE], ignore_attr = TRUE)
}

spotlight_fixture <- function() {
  fixture <- reference_mixture()
  genes <- paste0("Gene", c(1:40, 301:340, 601:640, 901:930))
  fixture$sce <- fixture$sce[genes, c(1:10, 41:50, 81:90)]
  fixture$spe <- fixture$spe[genes, ]
  fixture
}

spotlight_markers <- function() {
  data.frame(gene = paste0("Gene", c(1:40, 301:340, 601:640)),
    cluster = rep(c("A", "B", "C"), each = 40), mean.AUC = .95)
}
