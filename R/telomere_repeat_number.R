#' Count the number of telomeric sequences repeated in a genome sequence
#'
#' This function counts the occurrences of specified telomere repeat sequences (start and end) within a genome sequence,
#' splitting the sequence into windows of a specified size (default is 1 Mb). It returns a data table with the counts
#' of telomere repeats for each window of each chromosome in the input genome.
#'
#' @param fasta A DNA sequences read in with the readDNAStringSet function.
#' @param window A numeric value specifying the size of the window (in base pairs) for splitting the genome sequences.
#'        Default is 1e6 (1 million base pairs).
#' @param tel_start A character string representing the start sequence of the telomere. Default is "CCCTAAA".
#' @param tel_end A character string representing the end sequence of the telomere. Default is "TTTAGGG".
#'
#' @return A data table with the following columns:
#' \describe{
#'   \item{Chromosome}{The name of the chromosome.}
#'   \item{Window}{The window number (corresponding to the 1 Mb segments of the genome).}
#'   \item{CCCTAAA_Counts}{The count of the start telomere repeat sequence ("CCCTAAA") in the window.}
#'   \item{TTTAGGG_Counts}{The count of the end telomere repeat sequence ("TTTAGGG") in the window.}
#' }
#'
#' @examples
#' # Example usage with a list of genome sequences
#' fasta <- list(Chr1 = "ATGCATGCATGCATGCATGCATGC", Chr2 = "GGGCCCTAAAGGGCCCTAAA")
#' telomere_repeat_number(fasta, window = 4, tel_start = "CCCTAAA", tel_end = "TTTAGGG")
#'
#' @seealso \code{\link{create_windows}} and \code{\link{count_sequence}} for the functions used to split sequences
#' and count the telomere repeats, respectively.
#'
#' @importFrom data.table data.table rbindlist
#' @importFrom pbapply pblapply
#'
#' @export
telomere_repeat_number <- function(fasta, window = 1e6, tel_start = "CCCTAAA", tel_end = "TTTAGGG") {
  # Set telomeres
  telomere_start <- strrep(tel_start, 3)
  telomere_end <- strrep(tel_end, 3)

  # Use pblapply to loop through each element in 'fasta'
  results_list <- pblapply(seq_along(fasta), function(i) {
    seq_name <- names(fasta[i])
    sequence <- as.character(fasta[i])

    # Split the sequence into 1 Mb windows
    windows <- create_windows(genome = sequence, window_size = window)

    # Count occurrences of the first telomere sequence in each window
    telomere_start_counts <- sapply(windows, count_sequence, sequence = telomere_start)

    # Count occurrences of the last telomere sequence in each window
    telomere_end_counts <- sapply(windows, count_sequence, sequence = telomere_end)

    # Create a data table for the window counts
    dt_chromosome <- data.table(
      Chromosome = rep(seq_name, length(windows)),
      Window = seq_along(windows),
      CCCTAAA_Counts = telomere_start_counts,
      TTTAGGG_Counts = telomere_end_counts
    )

    return(dt_chromosome)
  })

  # Combine the list of data tables into one data table
  dt <- rbindlist(results_list)

  return(dt)
}
