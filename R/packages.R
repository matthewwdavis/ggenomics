# This function runs automatically when the package is loaded with `library(ggideo)`.
# It loads the required libraries so they are available to users.

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Loading ggideo and required packages...")

  # Load each package
  suppressMessages({
    library(data.table)
    library(tidyverse)
    library(scales)
    library(Biostrings)
    library(pbapply)
  })
}
