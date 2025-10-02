#' Manage python dependencies
#' according to: https://rstudio.github.io/reticulate/articles/python_dependencies.html#manual-configuration
#'
#' @name spacedeconvstartup
NULL

.onLoad <- function(libname, pkgname) {

  cli::cli_alert("checking spacedeconv environment and dependencies")

  envname <- getOption("omnideconv.conda_env", default = "r-omnideconv")

  cli::cli_alert("Using conda environment '{envname}'")

  if (!dir.exists(reticulate::miniconda_path())) {
    message("Setting python version in miniconda to be 3.10")
    Sys.setenv(RETICULATE_MINICONDA_PYTHON_VERSION = 3.10)
    message("Installing miniconda..")
    suppressMessages(reticulate::install_miniconda())
  }

  if (!(envname %in% reticulate::conda_list()$name)) {
    stop(cli::format_error("Conda environment '{envname}' not found."), call. = FALSE)
  }

  reticulate::use_miniconda(condaenv = envname, required = TRUE)

  # bug fix
  Csparse_validate <- "CsparseMatrix_validate"
}
