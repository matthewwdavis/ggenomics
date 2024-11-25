#' Create a Base ggplot Object for Genomic Data
#'
#' This function initializes a ggplot2 object for genomic data, allowing for the seamless integration of additional layers and customization. 
#' It is designed as a lightweight wrapper around [ggplot2::ggplot()] with additional flexibility for genomic data visualization.
#'
#' @param data A data frame or tibble containing the genomic data to be plotted.
#' @param mapping Default list of aesthetic mappings to be used for the plot. Must be created using [ggplot2::aes()] or [ggplot2::aes_string()]. Default is `NULL`.
#' @param ... Other arguments passed on to [ggplot2::ggplot()], such as `environment`.
#'
#' @return A ggplot object that can be further customized with additional layers, themes, and scales.
#'
#' @details
#' The `ggnom` function serves as a starting point for creating genomic plots. By wrapping the `ggplot` function, it simplifies workflows in genomic visualization 
#' while maintaining full compatibility with ggplot2 functionality.
#'
#' @examples
#' library(ggplot2)
#'
#' # Example dataset
#' df <- data.frame(
#'   Chromosome = factor(1:3),
#'   Position = c(100, 200, 300),
#'   Value = c(0.5, 0.7, 0.2)
#' )
#'
#' # Create a base genomic plot
#' ggnom(df, aes(x = Position, y = Value, color = Chromosome)) +
#'   geom_point() +
#'   labs(title = "Genomic Data Plot", x = "Position", y = "Value")
#'
#' @seealso [ggplot2::ggplot()], [ggplot2::aes()]
#'
#' @import ggplot2
#' @export
ggnom <- function(data, mapping = NULL, ...) {
  ggplot(data = data, mapping = mapping, ...)
}


