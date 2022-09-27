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
      reticulate::use_miniconda(condaenv = "r-reticulate", required = TRUE)
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
anndata_checkload <- function() {
  if (!python_available()) {
    base::message("Setting up python environment..")
    init_python()
    if (!python_available()) {
      base::stop(
        "Could not initiate miniconda python environment. Please set up manually with ",
        "init_python(python=your/python/version)"
      )
    }
  }
  if (!reticulate::py_module_available("anndata")) {
    anndata::install_anndata()
  }
}
