#' Create a Table for Telomere Plotting
#'
#' This function generates a table of chromosomes, their lengths, and the associated telomere data for plotting purposes.
#' The table includes genomic information, telomere repeat counts, and proper formatting to ensure correct chromosome ordering and visualization.
#'
#' @param genome A fasta file read in by `DNAStringSet` to be analyzed.
#' @param chr_names A character string indicating the prefix for chromosome names (e.g., `"Chr"`). Default is `"Chr"`.
#' @param string_remove A character string to remove from chromosome names. Default is `"_RagTag"`.
#' @param tel_start_seq A character string representing the telomere start sequence. Default is `"CCCTAAA"`.
#' @param tel_end_seq A character string representing the telomere end sequence. Default is `"TTTAGGG"`.
#' @param size_windows A numeric value specifying the window size (in base pairs) for telomere repeat counting. Default is `1e6`.
#' @param min_tel_count A numeric value specifying the minimum telomere repeat count to include in the final table. Default is `25`.
#' @param sample_name A character string representing the sample name. Default is `NULL`.
#'
#' @return A data frame with columns for chromosome names, chromosome lengths, and telomere repeat counts. Chromosome names are properly formatted and ordered.
#'
#' @details 
#' This function performs several tasks:
#' - Filters chromosome data based on the provided chromosome name prefix.
#' - Counts telomeric repeats in windows of the genome using specified start and end telomere sequences.
#' - Filters out chromosomes with fewer telomere repeats than the specified threshold (`min_tel_count`).
#' - Formats chromosome names (removes specific strings, removes leading zeros, and ensures correct ordering).
#' 
#' The resulting table can be used for plotting telomere data and chromosome lengths.
#'
#' @examples
#' \dontrun{
#' library(Biostrings)
#' genome <- readDNAStringSet("path/to/genome.fasta")
#' telomere_plotting_table(genome, chr_names = "Chr", string_remove = "_RagTag")
#' }
#'
#' @seealso [select_chr()], [telomere_repeat_number()], [genome_table()], [remove_string_chr()], [remove_lead_0s()], [remove_trailing()]
#'
#' @import data.table dplyr
#' @export
telomere_plotting_table <- function(genome, chr_names = "Chr", string_remove = "_RagTag", tel_start_seq = "CCCTAAA", tel_end_seq = "TTTAGGG",
                                    size_windows = 1e6, min_tel_count = 25, sample_name = NULL) {
  # Create table of contigs, chromosomes, and lengths
  length.table <- data.table(Chromosome = names(genome), Length = width(genome))
  
  # Filter for only chromosomes, based on starting string
  length.table <- select_chr(length.table, chr_string = chr_names)
  
  # Extract the size of the genome
  genome.size <- sum(length.table$Length)
  
  # Count telomeric sequence repeat
  tel_count.table <- telomere_repeat_number(fasta = genome, window = size_windows, tel_start = tel_start_seq, tel_end = tel_end_seq)
  
  # Filter to maintain telomeric counts over a certain threshold
  tel.table <- tel_count.table %>%
    filter(Forward_Counts >= min_tel_count | Reverse_Counts >= min_tel_count)
  
  # Create the larger table necessary for plotting
  plotting.table <- genome_table(length.table, tel.table, name = sample_name, genome_size = genome.size)
  
  # Remove string from names
  plotting.table$Chromosome <- remove_string_chr(plotting.table, remove_string = string_remove)
  
  # Remove leading 0s for proper ordering and plotting
  plotting.table <- remove_lead_0s(plotting.table, chr_string = chr_names)
  
  # Remove trailing strings. These is defined by a space
  plotting.table <- remove_trailing(plotting.table)
  
  # Set levels so that chromosomes are plotted in the proper order by number
  plotting.table$Chromosome <- factor(plotting.table$Chromosome,
                                      levels = unique(plotting.table$Chromosome)[
                                        order(sapply(unique(plotting.table$Chromosome), function(x) {
                                          # If the chromosome value is purely numeric, keep it as is
                                          if (grepl("^\\d+$", x)) {
                                            return(as.numeric(x))  # Convert numeric strings to numbers
                                          } else {
                                            # For non-numeric chromosomes, remove everything after the space
                                            return(as.character(gsub(" .*", "", x)))
                                          }
                                        }))
                                      ])
  
  return(plotting.table)
}
