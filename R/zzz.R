#' Manage python dependencies
#' according to: https://rstudio.github.io/reticulate/articles/python_dependencies.html#manual-configuration
#'
#' @name spacedeconvstartup
NULL

.onLoad <- function(libname, pkgname) {
  tmp <- utils::capture.output({
    # Make sure miniconda is installed
    if (!dir.exists(reticulate::miniconda_path())) {
      message("Setting python version in miniconda to be 3.8")
      Sys.setenv(RETICULATE_MINICONDA_PYTHON_VERSION = 3.8)
      message("Setting up miniconda environment..")
      suppressMessages(reticulate::install_miniconda())
    }


    # We ensure to have the r-spacedeconv env
    # if (!file.exists(reticulate::conda_python("r-reticulate"))) {
    if (!("r-omnideconv" %in% reticulate::conda_list()$name)) {
      reticulate::conda_create(envname = "r-omnideconv")
    }

    paths <- reticulate::conda_list()
    path <- paths[paths$name == "r-omnideconv", 2]
    if (.Platform$OS.type == "windows") {
      path <- gsub("\\\\", "/", path)
    }
    path.bin <- gsub("/envs/r-omnideconv/python.exe", "/library/bin", path)
    Sys.setenv(PATH = paste(path.bin, Sys.getenv()["PATH"], sep = ";"))
    Sys.setenv(RETICULATE_PYTHON = path)


    reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)
    reticulate::py_config()
    reticulate::configure_environment(pkgname, force = TRUE)

    if (!reticulate::py_module_available("anndata")) {
      anndata::install_anndata()
    }

  # if (!reticulate::py_module_available("python.app")) { # only MAC!
  #   reticulate::py_install("python.app", pip=TRUE)
  # }

    if (!reticulate::py_module_available("sklearn")) {
      reticulate::py_install("scikit-learn")
    }

    if (!reticulate::py_module_available("community")) {
      reticulate::py_install("python-louvain", pip = TRUE)
    }
    if (!reticulate::py_module_available("cell2location")) {
      reticulate::py_install("cell2location", pip = TRUE)
    }
  })
}
