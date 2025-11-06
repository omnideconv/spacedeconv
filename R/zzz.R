#' Manage python dependencies
#' according to: https://rstudio.github.io/reticulate/articles/python_dependencies.html#manual-configuration
#'
#' @name spacedeconvstartup
NULL

.onLoad <- function(libname, pkgname) {

  cli::cli_alert("Setting up spacedeconv environment..")

  Sys.setenv(RETICULATE_AUTOCONFIGURE = "FALSE")
  options(Seurat.object.assay.version = "v3")

  envname <- getOption("omnideconv.conda_env", default = "r-omnideconv")

  cli::cli_alert("Using conda environment '{envname}'")

  if (!dir.exists(reticulate::miniconda_path())) {
    cli::cli_alert("Setting python version in miniconda to be 3.10")
    Sys.setenv(RETICULATE_MINICONDA_PYTHON_VERSION = "3.10")
    cli::cli_alert("Installing miniconda..")
    reticulate::install_miniconda()
    cli::cli_alert("Miniconda installation complete")
  }

  if (!(envname %in% reticulate::conda_list()$name)) {
    stop(cli::format_error("Conda environment '{envname}' not found."), call. = FALSE)
  }

  reticulate::use_miniconda(condaenv = envname, required = TRUE)

  # bug fix
  Csparse_validate <- "CsparseMatrix_validate"
}
