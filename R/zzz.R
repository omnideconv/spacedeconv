#' Manage python dependencies
#' according to: https://rstudio.github.io/reticulate/articles/python_dependencies.html#manual-configuration
#'
#' @name spacedeconvstartup
NULL


.required_python_modules <- list(
  list(pypi = "python-igraph",   import = "igraph"),
  list(pypi = "leidenalg",       import = "leidenalg"),
  list(pypi = "python-louvain",  import = "community"),
  list(pypi = "scikit-learn",    import = "sklearn"),
  list(pypi = "scanpy",          import = "scanpy"),
  list(pypi = "cell2location",   import = "cell2location"),
  list(pypi = "scaden",          import = "scaden")
)


.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "This package uses Python via the 'reticulate' package.\n",
    "To install Python dependencies in your reticulate environment, run:\n",
    "    spacedeconv::setup_python_environment()\n"
  )

  if (!reticulate::py_available(initialize = TRUE)) {
    cli::cli_alert_warning("Python not yet initialized; skipping module check.")
    return()
  }

  py_ver <- as.character(reticulate::py_config()$version)
  if (utils::compareVersion(py_ver, "3.10") < 0) {
    cli::cli_alert_danger("Python >= 3.10 is required, but version {py_ver} was found.")
    return()
  }

  for (mod in .required_python_modules) {
    if (!check_python_module(mod$import)) {
      cli::cli_alert_danger("Python module missing: {mod$import}")
    } else {
      cli::cli_alert_success("Python module available: {mod$import}")
    }
  }
}


#' @export
setup_python_environment <- function() {
  pkgs <- unique(vapply(.required_python_modules, function(x) x$pypi, character(1)))
  reticulate::py_require(
    packages = pkgs,
    python_version = ">=3.10"
  )
}


check_python_module <- function(module) {
  code <- sprintf("
import importlib.util
result = importlib.util.find_spec('%s') is not None
", module)
  reticulate::py_run_string(code)
  reticulate::py$result
}