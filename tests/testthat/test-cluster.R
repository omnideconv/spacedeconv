data("spatial_data_3")

# Test if 'cluster' returns a SpatialExperiment object
# test_that("cluster returns a SpatialExperiment object", {
#   result <- cluster(spe = spacedeconv::preprocess(spatial_data_3, remove_mito = TRUE), spmethod = "expression")
#   expect_s4_class(result, "SpatialExperiment")
# })

# Test if 'cluster' handles null 'spe' input correctly
test_that("cluster handles null spe input correctly", {
  expect_error(cluster(spe = NULL))
})




spe <- readRDS(system.file("testdata", "spe.rds", package = "spacedeconv"))
spe <- spacedeconv::preprocess(spe, remove_mito = TRUE)
spe <- deconvolute(spe, method = "estimate")


# Test if 'cluster' adds clustering results to the SpatialExperiment object
test_that("cluster adds clustering results", {
  result <- cluster(spe = spe, method = "kmeans", spmethod = "estimate", nclusters = 3)
  expect_true("cluster_estimate_nclusters_3" %in% colnames(colData(result)))
})


test_that("get_cluster_features returns expected top features", {
  cluster <- cluster(spe = spe, method = "hclust", spmethod = "estimate", nclusters = 3)
  result <- get_cluster_features(cluster, clusterid = "cluster_estimate_nclusters_3", topn = 2, spmethod = "estimate")
  expect_equal(length(result[[1]]), 2)
})

test_that("get_cluster_features handles null spe input correctly", {
  expect_error(get_cluster_features(spe = NULL, clusterid = "clusterid"))
})

test_that("get_cluster_features handles null clusterid input correctly", {
  expect_error(get_cluster_features(spe, clusterid = NULL))
})
