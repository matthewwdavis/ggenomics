#' Telomere Plot Components for ggplot2
#'
#' This function creates a list of ggplot2 components for visualizing telomere size and chromosome length in a plot.
#'
#' @param chr_color A character string specifying the color of the chromosome segments. Default is `"#F8766D"`.
#' @param chr_size A numeric value specifying the thickness of the chromosome segments. Default is `6`.
#' @param tel_color A character string specifying the color of the telomere points. Default is `"black"`.
#' @param tel_shape An integer specifying the shape of the telomere points. Default is `16`.
#' @param legend_title A character string specifying the title for the legend representing telomere size. Default is `"Telomere Size (bp)"`.
#' @param text_size A numeric value specifying the base size for plot text, such as axis labels and titles. Default is `6`.
#' @param plot_title A character string specifying the title of the plot. Default is `NULL` (no title).
#' @param x_axis_title A character string specifying the title of the x-axis. Default is `NULL` (no title).
#' @param y_axis_title A character string specifying the title of the y-axis. Default is `"Chromosome Length"`.
#'
#' @return A list of ggplot2 components for building a telomere plot in ideogram style.
#'
#' @details 
#' The function returns a list of ggplot2 layers and components:
#' - Chromosome segments (`geom_segment`) are plotted from the telomere start positions to chromosome lengths.
#' - Telomere points (`geom_point`) are added at the start and end positions of the telomeres, with sizes proportional to telomere lengths.
#' - A customizable legend for telomere size.
#' - Optional titles for the plot and axes.
#' 
#' Users can add this list of components to a ggplot2 object to generate a telomere plot.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # Example dataset
#' df <- data.frame(
#'   Chromosome = factor(1:3),
#'   begin_telo_start = c(0, 0, 0),
#'   begin_telo_end = c(100, 200, 300),
#'   begin_telo_bp = c(100, 200, 300),
#'   end_telo_end = c(1500, 2000, 2500),
#'   end_telo_bp = c(500, 400, 300),
#'   Length = c(2000, 2500, 3000)
#' )
#'
#' # Telomere plot
#' ggplot(df) +
#'   geom_telplot(
#'     chr_color = "blue",
#'     tel_color = "red",
#'     plot_title = "Telomere Visualization",
#'     x_axis_title = "Chromosomes",
#'     y_axis_title = "Chromosome Length (bp)"
#'   )
#' }
#'
#' @seealso [ggplot2::geom_segment()], [ggplot2::geom_point()], [ggplot2::labs()], [ggplot2::theme()]
#'
#' @import ggplot2
#' @export
geom_telplot <- function(chr_color = "#F8766D", chr_size = 6,
                         tel_color = "black", tel_shape = 16, 
                         legend_title = "Telomere Size (bp)",
                         text_size = 6, plot_title = NULL, x_axis_title = NULL,
                         y_axis_title = "Chromosome Length") {
  list(
    geom_segment(aes(xend = Chromosome, y = begin_telo_start, yend = Length),
                 color = chr_color,
                 linewidth = chr_size,
                 lineend = "round"),
    geom_point(aes(x = Chromosome, y = begin_telo_end, size = ifelse(begin_telo_bp == 0, NA, begin_telo_bp)),
               shape = tel_shape,
               color = tel_color,
               na.rm = TRUE),
    geom_point(aes(x = Chromosome, y = end_telo_end, size = ifelse(end_telo_bp == 0, NA, end_telo_bp)),
               shape = tel_shape,
               color = tel_color,
               na.rm = TRUE),
    labs(y = y_axis_title, x = x_axis_title, size = legend_title, title = plot_title),
    theme_classic(base_size = text_size) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )
}

