% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ideogram_diploid.R
\name{ideogram_diploid}
\alias{ideogram_diploid}
\title{Create an ideogram plot for haplotype phased diploid genomic assemblies}
\usage{
ideogram_diploid(
  genome.table,
  plot_title = NULL,
  x_axis_title = NULL,
  y_axis_title = "Chromosome Length",
  legend_title = "Telomere Size (bp)",
  hap1_name = "Haplotype 1",
  hap2_name = "Haplotype 2",
  hap1_color = "dodgerblue2",
  hap2_color = "orangered",
  chr_size = 6,
  chr_distance = 0.7,
  tel_color = "black",
  tel_shape = 16,
  y_scale = 1e-06,
  y_scale_suffix = "Mb",
  legend_chr_size = 3,
  legend_pos = "right",
  legend_size = 0.25,
  text_size = 6
)
}
\arguments{
\item{genome.table}{A data frame containing genomic data, with columns for chromosome names, chromosome lengths, haplotypes,
and the start and end positions of telomere regions. (\code{begin_telo_start}, \code{begin_telo_end}, \code{end_telo_start},
\code{end_telo_end}, \code{begin_telo_bp}, and \code{end_telo_bp}).}

\item{plot_title}{A character string specifying the plot title. Default is \code{NULL}, in which case no title is added.}

\item{x_axis_title}{A character string specifying the x axis title. Default is \code{NULL}, in which case no title is added.}

\item{y_axis_title}{A character string specifying the y axis title. Default is "Chromosome Length".}

\item{legend_title}{A character string specifying the legend title. Default is "Telomere Presence".}

\item{hap1_name}{A character string specifying the name of the haplotype. Default is "Haplotype 1".}

\item{hap2_name}{A character string specifying the name of the haplotype. Default is "Haplotype 2".}

\item{hap1_color}{A character string specifying the color of the chromosome segments. Default is "dodgerblue2".}

\item{hap2_color}{A character string specifying the color of the chromosome segments. Default is "orangered".}

\item{chr_size}{A numeric value specifying the size (linewidth) of the chromosome segments. Default is 6.}

\item{chr_distance}{A numeric value specifying the distance between chromosome segments. Default is 0.7.}

\item{tel_color}{A character string specifying the color of the telomere points. Default is "black".}

\item{tel_shape}{A numeric value specifying the shape of the telomere points. Default is 16 (filled circle).}

\item{y_scale}{A numeric value for scaling the y-axis. Default is \code{1e-6} (for scaling the length to Mb).}

\item{y_scale_suffix}{A character string to append to the y-axis labels, usually a unit suffix like "Mb". Default is "Mb".}

\item{legend_chr_size}{An integer to change the size of the chromosomes in the legend. Default is 3.}

\item{legend_pos}{A character string or coordinates to specify the position of the legend. Default is "right".}

\item{legend_size}{A numeric value specifying the size of the legend keys. Default is 0.25.}

\item{text_size}{A numeric value specifying the size of text for the entire figure. Adjusts \code{base_size} ggplot function. Default is 6.}
}
\value{
A \code{ggplot} object representing the ideogram plot with chromosome and telomere information visualized.
}
\description{
This function generates a basic ideogram plot of chromosomes, including the chromosome length and the telomere regions.
It visualizes the telomeres by marking the start and end points, and adjusting the plot appearance according to
specified parameters such as color, size, and scaling. This function was built to use datasets generated with genome_table
}
\examples{
# Example usage
genome.table <- data.frame(
  Chromosome = c("Chr1", "Chr2"),
  Length = c(100000000, 200000000),
  begin_telo_start = c(0, 0),
  begin_telo_end = c(5000000, 10000000),
  end_telo_start = c(95000000, 190000000),
  end_telo_end = c(100000000, 200000000),
  begin_telo_bp = c(5000000, 10000000),
  end_telo_bp = c(5000000, 10000000)
)
ideogram_diploid(genome.table)

}
