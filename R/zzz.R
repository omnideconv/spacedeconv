#' Manage python dependencies
#' according to: https://rstudio.github.io/reticulate/articles/python_dependencies.html#manual-configuration
#'
#' @name spacedeconvstartup
NULL

.__conda_env_default_name__ <- "r-omnideconv_test"

.__required_python_version__ <- "3.10"

.required_python_modules <- list(
  list(pypi = "python-igraph",   import = "igraph",        install = "conda"),
  list(pypi = "leidenalg",       import = "leidenalg",     install = "conda"),
  list(pypi = "python-louvain",  import = "community",     install = "conda"),
  list(pypi = "scikit-learn",    import = "sklearn",       install = "conda"),
  list(pypi = "scanpy",          import = "scanpy",        install = "conda"),
  list(pypi = "cell2location",   import = "cell2location", install = "pip"),
  list(pypi = "scaden",          import = "scaden",        install = "pip")
)

.onAttach <- function(libname, pkgname) {

  check_reticulate_env()

  if(!all_py_modules_available()){
    packageStartupMessage(
      "This package uses Python via the 'reticulate' package.\n",
      "To setup (and install if needed) Python dependencies in your reticulate environment, run:\n",
      "    spacedeconv::setup_python_environment()\n",
      "Choose the right reticulate environment if you have already done this step."
    )
    display_py_modules_status()
  }

}

#' @export
setup_python_environment <- function() {
  
  reticulate::install_miniconda()

  envname <- .__conda_env_default_name__
  if (!reticulate::py_available(initialize = FALSE)) {
    if(!(envname %in% reticulate::conda_list()$name)){
      reticulate::conda_create(envname, python_version = .__required_python_version__)
    }
    reticulate::use_condaenv(envname, required = TRUE)
  }

  check_reticulate_env()

  # Separate packages by install method
  conda_pkgs <- vapply(
    .required_python_modules, 
    function(x) if (x$install == "conda") x$pypi else NULL, 
    character(1)
  )
  pip_pkgs <- vapply(
    .required_python_modules, 
    function(x) if (x$install == "pip") x$pypi else NULL, 
    character(1)
  )

  conda_pkgs <- conda_pkgs[nzchar(conda_pkgs)]
  pip_pkgs <- pip_pkgs[nzchar(pip_pkgs)]

  # Install conda packages
  if (length(conda_pkgs) > 0) {
    reticulate::conda_install(
      envname = NULL,
      packages = conda_pkgs,
      channel = "conda-forge"
    )
  }

  # Install pip packages
  if (length(pip_pkgs) > 0) {
    reticulate::py_install(
      packages = pip_pkgs,
      envname = NULL,
      pip = TRUE
    )
  }

  reticulate::py_config()
  display_py_modules_status()
}

check_python_module <- function(module) {
  code <- sprintf("
import importlib.util
result = importlib.util.find_spec('%s') is not None
", module)
  reticulate::py_run_string(code)
  reticulate::py$result
}

check_reticulate_env <- function() {
  if (reticulate::py_available(initialize = FALSE)) {
    cfg <- reticulate::py_config()
    py_ver <- as.character(cfg$version)
    envname <- if (is.null(cfg$condaenv)) "<unknown>" else cfg$condaenv

    if (is.null(cfg$conda)) {
      stop(sprintf(
        "Reticulate is using a Python environment that is NOT managed by Conda.\nThis package requires a conda environment for reticulate."
      ))
    }

    if (utils::compareVersion(py_ver, .__required_python_version__) != 0) {
      stop(sprintf(
        "python=%s is required, but reticulate was already initialized with Python version %s (conda env: %s).",
        .__required_python_version__, py_ver, envname
      ))
    }
  }
}

display_py_modules_status <- function() {
  if (reticulate::py_available(initialize = FALSE)) {
    for (mod in .required_python_modules) {
      if (!check_python_module(mod$import)) {
        cli::cli_alert_danger("Python module missing: {mod$import}")
      } else {
        cli::cli_alert_success("Python module available: {mod$import}")
      }
    }  
  } else {
    cli::cli_alert_warning("Python not yet initialized; skipping module check.")
  }
}

all_py_modules_available <- function() {
  if (!reticulate::py_available(initialize = FALSE)) {
    return(FALSE)
  }

  all_available <- TRUE

  for (mod in .required_python_modules) {
    if (!check_python_module(mod$import)) {
      all_available <- FALSE
      break
    }
  }

  return(all_available)
}