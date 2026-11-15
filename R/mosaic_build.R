#' Rebuild a mosaic Figure from Its Plots
#'
#' This function runs `mosaic build` on a mosaic document: every linked plot is read again from its
#' file, keeping its position and scale, and the figure is exported. Call it at the end of a script
#' that re-makes the plots (e.g. with [save_for_mosaic()]) so the finished figure follows the data
#' without opening mosaic. The document itself is left as it is unless `save = TRUE`.
#'
#' @param document The path of the mosaic document (.mosaic).
#' @param formats The formats to export: any of "pdf", "svg", "png", "tiff", "jpeg" and "eps". Default is "pdf".
#' @param output The path of the exported files, without an extension. Default is NULL, which names
#' them after the document.
#' @param dpi The resolution of exported images. Default is 600.
#' @param pages The numbers of the pages to export. Default is NULL, which exports every page.
#' @param journal The name of a journal whose panel-label style the export should use (e.g. "Nature").
#' Default is NULL, which keeps the labels as they are.
#' @param save If TRUE, the document is saved with the new plots too (the file as it was is kept
#' under File > Browse Versions in mosaic). Default is FALSE.
#' @param mosaic The mosaic command to run. Default is NULL, which looks for it with [find_mosaic()].
#' @param quiet If TRUE, mosaic's notes are not printed. Default is FALSE.
#'
#' @return The paths of the exported files, invisibly.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' mosaic_build("./Figures/main_figures.mosaic", formats = c("pdf", "png"))
#'
#' # One page, in a journal's label style, saved under another name
#' mosaic_build("./Figures/main_figures.mosaic", pages = 2, journal = "Nature",
#'              output = "./Figures/submission/figure_2")
#' }
#'
#' @seealso [save_for_mosaic()], [find_mosaic()]
#' @export
mosaic_build <- function(document, formats = "pdf", output = NULL, dpi = 600, pages = NULL,
                         journal = NULL, save = FALSE, mosaic = NULL, quiet = FALSE) {
  # 1. Check the document and find mosaic
  if (!file.exists(document)) {
    stop("There is no mosaic document at ", document)
  }
  if (is.null(mosaic)) {
    mosaic <- find_mosaic()
  }

  # 2. Put the command's arguments together
  build.arguments <- c("build", shQuote(normalizePath(document)), "-f", paste(formats, collapse = ","),
                       "--dpi", as.character(dpi))
  if (!is.null(output)) {
    build.arguments <- c(build.arguments, "-o", shQuote(path.expand(output)))
  }
  if (!is.null(pages)) {
    build.arguments <- c(build.arguments, "--pages", paste(pages, collapse = ","))
  }
  if (!is.null(journal)) {
    build.arguments <- c(build.arguments, "--journal", shQuote(journal))
  }
  if (save) {
    build.arguments <- c(build.arguments, "--save")
  }

  # 3. Run mosaic: the files it wrote come back on standard output, its notes and warnings on standard error
  notes.file <- tempfile(fileext = ".txt")
  on.exit(unlink(notes.file), add = TRUE)
  exported.files <- suppressWarnings(system2(mosaic[1], c(mosaic[-1], build.arguments), stdout = TRUE,
                                             stderr = notes.file))
  status <- attr(exported.files, "status")
  notes <- readLines(notes.file, warn = FALSE)

  # 4. Stop if mosaic failed, with what it said
  if (!is.null(status) && status != 0) {
    stop("mosaic build failed:\n", paste(c(notes, exported.files), collapse = "\n"))
  }

  # 5. Pass on mosaic's notes and warnings
  if (!quiet && length(notes) > 0) {
    message(paste(notes, collapse = "\n"))
  }

  return(invisible(as.character(exported.files)))
}

#' Find the mosaic Command
#'
#' This function finds how to run mosaic on this computer. It looks, in order, at the option
#' `ggenomics.mosaic`, the environment variable `MOSAIC`, the `mosaic` command on the PATH, the
#' mosaic app (Mac), and the conda environment mosaic's installer makes.
#'
#' @return A character vector: the program to run, followed by any arguments that go before mosaic's own.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' find_mosaic()
#'
#' # Tell ggenomics where mosaic is, if it is somewhere unusual
#' options(ggenomics.mosaic = "/opt/tools/mosaic/bin/mosaic")
#' }
#'
#' @export
find_mosaic <- function() {
  # 1. A choice made by hand: the option, then the environment variable
  chosen <- getOption("ggenomics.mosaic", default = Sys.getenv("MOSAIC", unset = ""))
  if (nzchar(chosen[1])) {
    return(chosen)
  }

  # 2. The mosaic command on the PATH
  on.path <- Sys.which("mosaic")
  if (nzchar(on.path)) {
    return(unname(on.path))
  }

  # 3. The Mac app: its own Python runs mosaic as a module
  for (applications.folder in c("/Applications", path.expand("~/Applications"))) {
    app.python <- file.path(applications.folder, "mosaic.app", "Contents", "MacOS", "python")
    if (file.exists(app.python)) {
      return(c(app.python, "-P", "-m", "mosaic"))
    }
  }

  # 4. The conda environment the installer makes (Mac and Linux: bin/mosaic; Windows: Scripts/mosaic.exe)
  home <- path.expand("~")
  conda.folders <- c(file.path(home, c("miniforge3", "mambaforge", "miniconda3", "anaconda3")),
                     "/opt/anaconda3", "/opt/miniconda3", "/opt/miniforge3",
                     file.path(Sys.getenv("USERPROFILE", unset = home), "miniforge3"))
  for (conda.folder in conda.folders) {
    for (command in c(file.path("envs", "mosaic", "bin", "mosaic"),
                      file.path("envs", "mosaic", "Scripts", "mosaic.exe"))) {
      environment.command <- file.path(conda.folder, command)
      if (file.exists(environment.command)) {
        return(environment.command)
      }
    }
  }

  stop("mosaic was not found. Install it, or say where it is with options(ggenomics.mosaic = \"/path/to/mosaic\")")
}
