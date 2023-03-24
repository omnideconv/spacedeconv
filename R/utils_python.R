#' Initiates python environment
#'
#' @param python (optional) If own python should be used please indicate it's binaries
#'
init_python <- function(python = NULL) {
  if (!reticulate::py_available()) {
    if (is.null(python)) {
      if (!dir.exists(reticulate::miniconda_path())) {
        message("Setting python version in miniconda to be 3.8")
        Sys.setenv(RETICULATE_MINICONDA_PYTHON_VERSION = 3.8)
        message("Setting up miniconda environment..")
        suppressMessages(reticulate::install_miniconda())
      }
      check_env()
      reticulate::use_miniconda(condaenv = "spacedeconv", required = TRUE)
      config <- reticulate::py_config()
      if (!python_available()) {
        message("Python not available")
        print(config)
        message(
          "Please indicate your version of python calling init_python(python=your/python)"
        )
      }
    } else {
      reticulate::use_python(python = python)
      reticulate::py_config()
    }
  }
}

check_env <- function() {
  if (!"spacedeconv" %in% reticulate::conda_list()$name) {
    reticulate::conda_create(envname = "spacedeconv", packages = c("pip", "numpy"))
    cli::cli_alert_info("created conda environment 'spacedeconv'")
  }
}


#' Checks if python is available in environment
#'
#' @return boolean
#'
python_available <- function() {
  return(reticulate::py_available())
}

#' Checks if anndata package is loaded
#'
#' If called and python environment is not set up, this is realized. Else, it checks if the anndata
#' package is loaded, and if not, it does this.
#'
# anndata_checkload <- function() {
#   if (!python_available()) {
#     base::message("Setting up python environment..")
#     init_python()
#     if (!python_available()) {
#       base::stop(
#         "Could not initiate miniconda python environment. Please set up manually with ",
#         "init_python(python=your/python/version)"
#       )
#     }
#   }
#   if (!reticulate::py_module_available("anndata")) {
#     anndata::install_anndata()
#   }
# }

#' #' Checks if metacells in installed
#' metacells_checkload <- function() {
#'   if (!python_available()) {
#'     base::message("Setting up python environment..")
#'     init_python()
#'     if (!python_available()) {
#'       base::stop(
#'         "Could not initiate miniconda python environment. Please set up manually with ",
#'         "init_python(python=your/python/version)"
#'       )
#'     }
#'   }
#'   if (!reticulate::py_module_available("metacells")) {
#'     reticulate::py_install("metacells", pip = TRUE)
#'   }
#' }

#' #' Checks if cell2location in installed
#' cell2location_checkload <- function() {
#'   if (!python_available()) {
#'     base::message("Setting up python environment..")
#'     init_python()
#'     if (!python_available()) {
#'       base::stop(
#'         "Could not initiate miniconda python environment. Please set up manually with ",
#'         "init_python(python=your/python/version)"
#'       )
#'     }
#'   }
#'   if (!reticulate::py_module_available("cell2location")) {
#'     reticulate::py_install("scanpy", pip = TRUE)
#'     reticulate::py_install("cell2location", pip = TRUE)
#'   }
#' }


# install_giotto_python <- function() {
#   if (!python_available()) {
#     base::message("Setting up python environment..")
#     init_python()
#     if (!python_available()) {
#       base::stop(
#         "Could not initiate miniconda python environment. Please set up manually with ",
#         "init_python(python=your/python/version)"
#       )
#     }
#   }
#   reticulate::py_install(c("python-igraph", "leidenalg", "python-louvain", "scikit-learn"))
#   reticulate::py_install("python.app", pip = TRUE)
# }

#' Install all python packages
#'
#' This makes sure a valid python installation exists and all needed packages are pulled and installed
#'
#' @export
# install_all_python <- function() {
#   init_python()
#   anndata_checkload()
#   metacells_checkload()
#   install_giotto_python()
#   cell2location_checkload()
#   omnideconv::install_all_python()
# }
