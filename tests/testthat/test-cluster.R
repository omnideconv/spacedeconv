test_that("cluster rejects the obsolete data argument instead of ignoring it", {
  spe <- result_spe()
  for (value in c("expression", "progeny")) {
    expect_error(cluster(spe, data = value), "Use 'spmethod'", fixed = TRUE)
  }
  expect_error(cluster(spe, spmethod = "quantiseq", data = "progeny"),
               "Use 'spmethod'", fixed = TRUE)
})

test_that("both result clustering methods recover separated groups", {
  local_test_state()
  spe <- result_spe()
  for (method in c("kmeans", "hclust")) {
    out <- cluster(spe, method = method, spmethod = "quantiseq", nclusters = 2,
                   dist_method = "euclidean")
    labels <- out$cluster_quantiseq_nclusters_2
    expect_length(unique(labels[1:3]), 1)
    expect_length(unique(labels[4:6]), 1)
    expect_false(labels[1] == labels[4])
    expect_preserved(out, spe)
    features <- get_cluster_features(out, "cluster_quantiseq_nclusters_2",
      topn = 1, spmethod = "quantiseq", zscore = TRUE)
    expect_setequal(unlist(lapply(features, names)), c("quantiseq_A", "quantiseq_B"))
  }
  expect_error(cluster(NULL), "null")
  expect_error(cluster(spe, spmethod = "unknown"))
  expect_error(get_cluster_features(NULL, "cluster"), "null")
})

test_that("expression clustering returns labels for the original spots", {
  local_test_state()
  spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv", mustWork = TRUE))
  spe <- preprocess(spe, min_umi = 0, remove_mito = TRUE)
  out <- cluster(spe, spmethod = "expression", pca_dim = 1:5, clusres = .5)
  expect_length(out$cluster_expression_res_0.5, ncol(spe))
  expect_false(anyNA(out$cluster_expression_res_0.5))
  expect_preserved(out, spe)
})
