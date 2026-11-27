#!/usr/bin/env Rscript
# Run from the repository root. This runs exactly the same suite as devtools::test().
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1]] else "test-results"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output)
options(testthat.progress.max_fails = Inf)
# Avoid accidental use of all CPUs on a shared machine.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
started <- Sys.time()
results <- devtools::test(reporter = "summary", stop_on_failure = FALSE)
saveRDS(results, file.path(output, "results.rds"))
issues <- unlist(lapply(results, function(test) {
  unlist(lapply(test$results, function(result) {
    if (inherits(result, "expectation_success")) return(NULL)
    paste(test$test, class(result)[1], conditionMessage(result), sep = "\n")
  }))
}))
writeLines(if (is.null(issues)) character() else issues, file.path(output, "issues.txt"))
table <- as.data.frame(results)
columns <- intersect(c("file", "test", "nb", "failed", "error", "warning", "skipped", "passed", "real"), names(table))
utils::write.csv(table[, columns], file.path(output, "results.csv"), row.names = FALSE)
sink(file.path(output, "session-info.txt"))
print(sessionInfo())
if (reticulate::py_available(initialize = FALSE)) print(reticulate::py_config())
cat("\nElapsed seconds:", as.numeric(difftime(Sys.time(), started, units = "secs")), "\n")
sink()
cat("\nResults:", output, "\n")
print(colSums(table[, c("failed", "error", "warning", "skipped", "passed")]))
if (any(table$skipped)) cat("INCOMPLETE: skipped tests are not verified. See the test output for prerequisites.\n")
# Distinguish incomplete coverage from a fully verified run in automation.
status <- if (any(table$failed > 0 | table$error)) 1L else if (any(table$skipped)) 2L else 0L
quit(status = status)
