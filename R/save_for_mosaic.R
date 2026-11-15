#' Save a Plot for a mosaic Figure
#'
#' This function saves a plot as a PDF, a PNG and an SVG with the same name and size, the three copies
#' mosaic places as one picture. Each file is written under a temporary name and only takes its real
#' name once it is complete, so a mosaic window that has the figure open never reads a half-written
#' file. mosaic notices the new files within a couple of seconds and updates the placed picture where
#' it is, keeping its position and scale.
#'
#' The plot is saved exactly as given: set its theme, margins and size before calling this function.
#' Backgrounds are transparent in all three formats.
#'
#' @param plot A ggplot object (or anything else `ggplot2::ggsave()` can draw).
#' @param stem The path of the files to write, without an extension (e.g. "./Figures/depth/chr1_depth").
#' The folder is created if it does not exist.
#' @param width The width of the saved figure, in `units`.
#' @param height The height of the saved figure, in `units`.
#' @param units The units of `width` and `height`: "in" (default), "cm", "mm" or "pt".
#' @param dpi The resolution of the PNG. Default is 600.
#' @param whole_points If TRUE (default), the size is rounded up to whole points (1/72 in). R's `pdf()`
#' cuts a page down to whole points, so without this the PDF comes out slightly smaller than the PNG
#' and SVG, and mosaic does not treat them as copies of one picture.
#' @param formats The formats to save. Default is c("pdf", "png", "svg"). The SVG needs the svglite package.
#' @param quiet If TRUE, no message is printed. Default is FALSE.
#'
#' @return The paths of the files written, invisibly.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' depth.plot <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg)) +
#'   ggplot2::geom_point() +
#'   ggplot2::theme_classic(base_size = 6, base_family = "Helvetica")
#' save_for_mosaic(depth.plot, "./Figures/demo/weight_mpg", width = 60, height = 45, units = "mm")
#'
#' # Re-make every export of the figure the plot is placed in
#' mosaic_build("./Figures/main_figures.mosaic", formats = c("pdf", "png"))
#' }
#'
#' @seealso [mosaic_build()]
#' @export
save_for_mosaic <- function(plot, stem, width, height, units = "in", dpi = 600,
                            whole_points = TRUE, formats = c("pdf", "png", "svg"), quiet = FALSE) {
  # 1. Check the formats and the size
  formats <- tolower(formats)
  unknown.formats <- setdiff(formats, c("pdf", "png", "svg"))
  if (length(unknown.formats) > 0) {
    stop("save_for_mosaic() saves pdf, png and svg; not ", paste(unknown.formats, collapse = ", "))
  }
  if ("svg" %in% formats && !requireNamespace("svglite", quietly = TRUE)) {
    stop("Saving the SVG needs the svglite package: install.packages(\"svglite\")")
  }
  to_inches <- c("in" = 1, "cm" = 1 / 2.54, "mm" = 1 / 25.4, "pt" = 1 / 72)
  if (!units %in% names(to_inches)) {
    stop("units must be one of: ", paste(names(to_inches), collapse = ", "))
  }
  width.in <- width * to_inches[[units]]
  height.in <- height * to_inches[[units]]

  # 2. Round the size up to whole points, plus a hair so pdf() does not cut it a point short
  if (whole_points) {
    width.in <- (ceiling(width.in * 72 - 1e-6) + 0.001) / 72
    height.in <- (ceiling(height.in * 72 - 1e-6) + 0.001) / 72
  }

  # 3. Make the folder
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)

  # 4. Save each format under a temporary name in the same folder, then give it its real name in one step
  written.files <- character(0)
  for (format in formats) {
    final.file <- paste0(stem, ".", format)
    temporary.file <- paste0(stem, ".saving-", Sys.getpid(), ".", format)
    on.exit(unlink(temporary.file), add = TRUE)
    if (format == "png") {
      ggplot2::ggsave(temporary.file, plot, width = width.in, height = height.in, units = "in", dpi = dpi,
                      bg = "transparent")
    } else {
      ggplot2::ggsave(temporary.file, plot, width = width.in, height = height.in, units = "in",
                      bg = "transparent")
    }
    if (!file.rename(temporary.file, final.file)) {
      stop("Could not write ", final.file)
    }
    written.files <- c(written.files, final.file)
  }

  # 5. Say what was saved
  if (!quiet) {
    message(sprintf("Saved %s (%s), %s x %s pt", stem, paste0(".", formats, collapse = " "),
                    format(round(width.in * 72, 1)), format(round(height.in * 72, 1))))
  }

  return(invisible(written.files))
}
