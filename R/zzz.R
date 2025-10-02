#' Manage python dependencies
#' according to: https://rstudio.github.io/reticulate/articles/python_dependencies.html#manual-configuration
#'
#' @name spacedeconvstartup
NULL

.onLoad <- function(libname, pkgname) {
  cli::cli_alert("checking spacedeconv environment and dependencies")

  envname <- getOption("omnideconv.conda_env", default = "r-omnideconv")

  cli::cli_alert("Using conda environment '{envname}'")

  temp_file <- tempfile()
  sink(temp_file)

  invisible({
    suppressMessages({
      suppressWarnings({
        # tmp <- reticulate::py_capture_output({
        tmp2 <- utils::capture.output({
          # reticulate::py_capture_output({
          # Make sure miniconda is installed
          if (!dir.exists(reticulate::miniconda_path())) {
            message("Setting python version in miniconda to be 3.10")
            Sys.setenv(RETICULATE_MINICONDA_PYTHON_VERSION = 3.10)
            message("Setting up miniconda environment..")
            suppressMessages(reticulate::install_miniconda())
          }


          # We ensure to have the r-spacedeconv env
          # if (!file.exists(reticulate::conda_python("r-reticulate"))) {
          if (!(envname %in% reticulate::conda_list()$name)) {
            reticulate::conda_create(envname = envname)
          }

          paths <- reticulate::conda_list()
          path <- paths[paths$name == envname, 2]
          if (.Platform$OS.type == "windows") {
            path <- gsub("\\\\", "/", path)
          }
          path.bin <- gsub(paste0("/envs/", envname, "/python.exe"), "/library/bin", path)
          Sys.setenv(PATH = paste(path.bin, Sys.getenv()["PATH"], sep = ";"))
          Sys.setenv(RETICULATE_PYTHON = path)


          reticulate::use_miniconda(condaenv = envname, required = FALSE)
          reticulate::py_config()
          reticulate::configure_environment(pkgname, force = FALSE)
        })
      })
    })
  })



  # bug fix
  Csparse_validate <- "CsparseMatrix_validate"

  sink(NULL)
  unlink(temp_file)
}
