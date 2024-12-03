.onAttach <- function(libname, pkgname) {
  suppressPackageStartupMessages({
    library(dplyr, warn.conflicts = FALSE)
    library(flowCore, warn.conflicts = FALSE)
  })
  packageStartupMessage("ANPELA package loaded successfully.")
}
