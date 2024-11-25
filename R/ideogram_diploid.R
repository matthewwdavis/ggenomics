#' Create an ideogram plot for haplotype phased diploid genomic assemblies
#'
#' This function generates a basic ideogram plot of chromosomes, including the chromosome length and the telomere regions.
#' It visualizes the telomeres by marking the start and end points, and adjusting the plot appearance according to
#' specified parameters such as color, size, and scaling. This function was built to use datasets generated with genome_table
#'
#' @param genome.table A data frame containing genomic data, with columns for chromosome names, chromosome lengths, haplotypes, 
#'        and the start and end positions of telomere regions. (`begin_telo_start`, `begin_telo_end`, `end_telo_start`,
#'        `end_telo_end`, `begin_telo_bp`, and `end_telo_bp`).
#' @param plot_title A character string specifying the plot title. Default is `NULL`, in which case no title is added.
#' @param x_axis_title A character string specifying the x axis title. Default is `NULL`, in which case no title is added.
#' @param y_axis_title A character string specifying the y axis title. Default is "Chromosome Length".
#' @param legend_title A character string specifying the legend title. Default is "Telomere Presence".
#' @param hap1_color A character string specifying the color of the chromosome segments. Default is "dodgerblue2".
#' @param hap2_color A character string specifying the color of the chromosome segments. Default is "orangered".
#' @param hap1_name A character string specifying the name of the haplotype. Default is "Haplotype 1".
#' @param hap2_name A character string specifying the name of the haplotype. Default is "Haplotype 2".
#' @param chr_size A numeric value specifying the size (linewidth) of the chromosome segments. Default is 6.
#' @param chr_distance A numeric value specifying the distance between chromosome segments. Default is 0.7.
#' @param tel_color A character string specifying the color of the telomere points. Default is "black".
#' @param tel_shape A numeric value specifying the shape of the telomere points. Default is 16 (filled circle).
#' @param y_scale A numeric value for scaling the y-axis. Default is `1e-6` (for scaling the length to Mb).
#' @param y_scale_suffix A character string to append to the y-axis labels, usually a unit suffix like "Mb". Default is "Mb".
#' @param legend_chr_size An integer to change the size of the chromosomes in the legend. Default is 3.
#' @param legend_pos A character string or coordinates to specify the position of the legend. Default is "right".
#' @param legend_size A numeric value specifying the size of the legend keys. Default is 0.25.
#' @param text_size A numeric value specifying the size of text for the entire figure. Adjusts `base_size` ggplot function. Default is 6.
#'
#' @return A `ggplot` object representing the ideogram plot with chromosome and telomere information visualized.
#'
#' @examples
#' # Example usage
#' genome.table <- data.frame(
#'   Chromosome = c("Chr1", "Chr2"),
#'   Length = c(100000000, 200000000),
#'   begin_telo_start = c(0, 0),
#'   begin_telo_end = c(5000000, 10000000),
#'   end_telo_start = c(95000000, 190000000),
#'   end_telo_end = c(100000000, 200000000),
#'   begin_telo_bp = c(5000000, 10000000),
#'   end_telo_bp = c(5000000, 10000000)
#' )
#' ideogram_diploid(genome.table)
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_point labs theme_classic theme element_text scale_y_continuous
#' @importFrom scales label_number
#' @importFrom dplyr filter
#'
#' @export
ideogram_diploid <- function(genome.table, plot_title = NULL, x_axis_title = NULL, y_axis_title = "Chromosome Length",
                             legend_title = "Telomere Size (bp)", hap1_name = "Haplotype 1", hap2_name = "Haplotype 2",
                             hap1_color = "dodgerblue2", hap2_color = "orangered", chr_size = 6, chr_distance = 0.7,
                             tel_color = "black", tel_shape = 16, y_scale = 1e-6, y_scale_suffix = "Mb",
                             legend_chr_size = 3, legend_pos = "right", legend_size = 0.25, text_size = 6) {
  
  # Plot the two haplotypes
  p <- genome.table %>%
    ggplot(aes(x = as.factor(Chromosome), y = Length)) +
    geom_segment(aes(y = begin_telo_start, yend = Length, color = Hap),
                 position = position_dodge(width = chr_distance),  # Use calculated dodge width
                 size = chr_size, 
                 lineend = "round") +
    geom_point(aes(x = Chromosome, y = begin_telo_end,
                   size = ifelse(begin_telo_bp == 0, NA, begin_telo_bp),
                   fill = Hap,
                   na.rm = TRUE),
               shape = tel_shape, color = tel_color, position = position_dodge(width = chr_distance)) +
    geom_point(aes(x = Chromosome, y = end_telo_end,
                   size = ifelse(end_telo_bp == 0, NA, end_telo_bp),
                   fill = Hap,
                   na.rm = TRUE),
               shape = tel_shape, color = tel_color, position = position_dodge(width = chr_distance)) +
    scale_color_manual(name = "Haplotype",
                       values = c(hap1 = hap1_color, hap2 = hap2_color),
                       labels = c(hap1 = hap1_name, hap2 = hap2_name)) +
    scale_y_continuous(
      labels = label_number(scale = y_scale, suffix = y_scale_suffix)) +
    labs(y = y_axis_title, x = x_axis_title, size = legend_title) +
    guides(color = guide_legend(title = NULL, override.aes = list(linewidth = legend_chr_size)),
           size = guide_legend(title = legend_title),
           fill = "none") +
    theme_classic(base_size = text_size) +
    theme(legend.position = legend_pos,
          legend.key.size = unit(legend_size, "cm"),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}
