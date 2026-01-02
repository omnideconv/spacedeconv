#' Manage python dependencies
#' according to: https://rstudio.github.io/reticulate/articles/python_dependencies.html#manual-configuration
#'
#' @name spacedeconvstartup
NULL

.onLoad <- function(libname, pkgname) {
  cli::cli_alert("Setting up spacedeconv environment..")

  Sys.setenv(RETICULATE_AUTOCONFIGURE = "FALSE")
  options(Seurat.object.assay.version = "v3")

  if (!dir.exists(reticulate::miniconda_path())) {
    cli::cli_alert("Setting python version in miniconda to be 3.10")
    Sys.setenv(RETICULATE_MINICONDA_PYTHON_VERSION = "3.10")
    cli::cli_alert("Installing miniconda..")
    reticulate::install_miniconda()
    cli::cli_alert("Miniconda installation complete")
  }

  default_envname <- "spacedeconv-env"
  envname <- getOption("spacedeconv.conda_env", default = default_envname)

  conda_envs <- reticulate::conda_list()
  if (!(envname %in% conda_envs$name)) {
    stop(cli::format_error("Conda environment '{envname}' not found."), call. = FALSE)
  }

  cli::cli_alert("Using conda environment '{envname}'")

  reticulate::use_miniconda(condaenv = envname, required = TRUE)

  # bug fix
  Csparse_validate <- "CsparseMatrix_validate"
}
