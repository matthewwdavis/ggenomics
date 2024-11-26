#' Create a Base ggplot Object for Genomic Data
#'
#' This function initializes a ggplot2 object for genomic data, allowing for the seamless integration of additional layers and customization. 
#' It is designed as a lightweight wrapper around [ggplot2::ggplot()] with additional flexibility for genomic data visualization.
#'
#' @param data A data frame or tibble containing the genomic data to be plotted.
#' @param mapping Default list of aesthetic mappings to be used for the plot. Must be created using [ggplot2::aes()] or [ggplot2::aes_string()]. Default is `NULL`.
#' @param plot A character string specifying a predefined plot type. Options include:
#'   - `"telplot"`: Sets mapping to `aes(x = Chromosome, y = begin_telo_start, yend = Length)`.
#' @param ... Other arguments passed on to [ggplot2::ggplot()], such as `environment`.
#'
#' @return A ggplot object that can be further customized with additional layers, themes, and scales.
#'
#' @details
#' The `ggenom` function serves as a starting point for creating genomic plots. By wrapping the `ggplot` function, it simplifies workflows in genomic visualization 
#' while maintaining full compatibility with ggplot2 functionality.
#'
#' @examples
#' library(ggplot2)
#'
#' # Example dataset
#' df <- data.frame(
#'   Chromosome = factor(1:3),
#'   begin_telo_start = c(100, 200, 300),
#'   Length = c(1000, 1500, 1200)
#' )
#'
#' # Using the plot argument
#' ggenom(df, plot = "telplot") +
#'   geom_telplot()
#'
#' @seealso [ggplot2::ggplot()], [ggplot2::aes()]
#'
#' @import ggplot2
#' @export
ggenom <- function(data, mapping = NULL, plot = NULL, ...) {
  # Define mapping presets
  plot_mappings <- list(
    telplot = aes(x = Chromosome, y = begin_telo_start, yend = Length)
  )
  
  # Automatically set mapping if plot is specified
  if (!is.null(plot)) {
    if (plot %in% names(plot_mappings)) {
      mapping <- plot_mappings[[plot]]
    } else {
      stop("Invalid plot type. Available options are: ", paste(names(plot_mappings), collapse = ", "))
    }
  }
  
  # Create the ggplot object
  ggplot(data = data, mapping = mapping, ...)
}
