sce <- readRDS(system.file("testdata", "sce.rds", package = "spacedeconv"))
library(SingleCellExperiment)

# Test the subsetSCE function
testthat::test_that("subsetSCE subsets the SingleCellExperiment correctly", {
  # Call the subsetSCE function
  result <- subsetSCE(sce, cell_type_col = "celltype_major", scenario = "even", ncells = 500, notEnough = "all", seed = 12345)

  # Check if the resulting object is a SingleCellExperiment
  expect_true(inherits(result, "SingleCellExperiment"))

  # Check the number of cells in the resulting object
  expect_equal(ncol(result), 250)
})


# Test case when the input SingleCellExperiment is NULL
testthat::test_that("subsetSCE throws an error when SingleCellExperiment is missing", {
  # Call the subsetSCE function with NULL input
  expect_error(
    subsetSCE(NULL, cell_type_col = "celltype_major", scenario = "even", ncells = 100, notEnough = "asis", seed = 12345),
    "SingleCellExperiment missing but required"
  )
})

# Test case when row and column names are not set in the SingleCellExperiment
testthat::test_that("subsetSCE throws an error when row or column names are not set", {
  # Create a mock SingleCellExperiment object without row or column names
  sce <- SingleCellExperiment(assays = list(counts = matrix(rpois(5000, lambda = 10), nrow = 100, ncol = 50)))

  # Call the subsetSCE function with the mock SingleCellExperiment
  expect_error(
    subsetSCE(sce, cell_type_col = "celltype_major", scenario = "even", ncells = 100, notEnough = "asis", seed = 12345),
    "Rownames or colnames are not set for SingleCellExperiment but are required"
  )
})


###########################

# Test: Verify if 'object' parameter is null or missing
test_that("Check for null or missing 'object' parameter", {
  expect_error(normalize(NULL), "Parameter 'object' is null or missing")
})

# Test: Verify if rownames and colnames are set
test_that("Check for missing rownames or colnames", {
  # Create an 'sce' object with missing rownames or colnames
  sce_missing_names <- SingleCellExperiment(
    assays = list(counts = matrix(1:12, nrow = 3))
  )
  expect_error(normalize(sce_missing_names), "Rownames or colnames not set")
})

# Test: Verify normalization method "cpm"
test_that("Check normalization using 'cpm' method", {
  # Create an 'sce' object for testing
  normalized_sce <- normalize(sce, method = "cpm")

  # Check if the assay "cpm" is added
  expect_true(c("cpm") %in% assayNames(normalized_sce))

  # Check if the assay "cpm" is of class "dgCMatrix"
  expect_true(is(assay(normalized_sce, "cpm"), "dgCMatrix"))
})

# Test: Verify normalization method "logcpm"
test_that("Check normalization using 'logcpm' method", {
  # Create an 'sce' object for testing
  normalized_sce <- normalize(sce, method = "logcpm")

  # Check if the assay "logcpm" is added
  expect_true("logcpm" %in% assayNames(normalized_sce))

  # Check if the assay "logcpm" is of class "dgCMatrix"
  expect_true(is(assay(normalized_sce, "logcpm"), "dgCMatrix"))
})

# Test: Verify return value of the function
test_that("Check return value of normalize function", {
  # Create an 'sce' object for testing
  normalized_sce <- normalize(sce)

  # Check if the returned object is of class "SingleCellExperiment"
  expect_true(is(normalized_sce, "SingleCellExperiment"))
})



############################


# Test: Verify if 'object' parameter is null or missing
test_that("Check for null or missing 'object' parameter", {
  expect_error(preprocess(NULL), "Please provide an object")
})

# Test: Verify return value of the function
test_that("Check return value of preprocess function", {
  preprocessed_sce <- preprocess(sce, remove_mito = TRUE)

  # Check if the returned object is of class "SingleCellExperiment"
  expect_true(is(preprocessed_sce, "SingleCellExperiment"))
})
