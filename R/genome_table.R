#' Create a table of genomic features with calculations
#'
#' This function generates a table combining information from a FASTA length table and a telomere table,
#' while calculating additional telomere-related features (such as the size of the telomeres at both ends).
#' It also normalizes the total telomere size relative to the total genome size.
#'
#' @param fasta_lengths_table A data frame containing the FASTA length information. Columns must be "Chromosome" and "Length".
#' @param telomere_table A data frame containing telomere-related data, including "Chromosome" and counts of telomere repeats (e.g., "CCCTAAA_Counts" and "TTTAGGG_Counts"). Generated from function telomere_repeat_number
#' @param name A character string representing the name or identifier of the genome source.
#' @param genome_size A numeric value representing the total size of the genome.
#'
#' @return A data frame with additional columns related to telomere sizes and their normalization. The new columns include:
#' \describe{
#'   \item{begin_telo_bp}{The length of the telomere at the beginning of the chromosome, in base pairs.}
#'   \item{end_telo_bp}{The length of the telomere at the end of the chromosome, in base pairs.}
#'   \item{begin_telo_start}{The start position of the beginning telomere, which is always 0.}
#'   \item{begin_telo_end}{The end position of the beginning telomere, calculated based on the `begin_telo_bp`.}
#'   \item{end_telo_start}{The start position of the ending telomere, calculated based on the chromosome length and `end_telo_bp`.}
#'   \item{end_telo_end}{The end position of the ending telomere, which is equal to the chromosome length.}
#'   \item{total_telo_bp}{The combined length of both telomeres in base pairs.}
#'   \item{normalized_total_telo_size}{The total telomere size normalized by the genome size.}
#'   \item{Source}{The name or identifier of the genome source (as provided by `name`).}
#' }
#'
#' @examples
#' # Example usage
#' fasta_lengths_table <- data.frame(Chromosome = c("Chr1", "Chr2"), Length = c(1000000, 2000000))
#' telomere_table <- data.frame(Chromosome = c("Chr1", "Chr2"), CCCTAAA_Counts = c(200, 300), TTTAGGG_Counts = c(150, 250))
#' genome_size <- 3000000
#' name <- "MyGenome"
#' genome_table(fasta_lengths_table, telomere_table, name, genome_size)
#'
#' @importFrom dplyr left_join select mutate
#'
#' @export
genome_table <- function(fasta_lengths_table, telomere_table, name, genome_size){
  table <- left_join(fasta_lengths_table, telomere_table, by = "Chromosome") %>%
    select(!Window) %>%
    mutate(begin_telo_bp = CCCTAAA_Counts*21) %>%
    mutate(end_telo_bp = TTTAGGG_Counts*21) %>%
    mutate(begin_telo_start = 0) %>%
    mutate(begin_telo_end = begin_telo_bp) %>%
    mutate(end_telo_start = Length - end_telo_bp) %>%
    mutate(end_telo_end = Length) %>%
    mutate(total_telo_bp = begin_telo_bp + end_telo_bp) %>%
    mutate(normalized_total_telo_size = total_telo_bp/genome_size) %>%
    mutate(Source = name)

  return(table)
}
