#' Create windows of specified size from the fasta file of a genome sequence
#'
#' This function splits a genome sequence into windows of a specified size. It uses the `substring`
#' function to extract substrings from the genome at regular intervals defined by the `window_size` parameter.
#'
#' @param genome A character string representing the genome sequence to split into windows. Should be a fasta read in with readDNAStringSet
#' @param window_size A numeric value representing the size of each window. Default is 1e6 (1 million).
#'
#' @return A character vector where each element is a window (substring) of the genome sequence.
#'
#' @examples
#' # Example usage
#' genome <- "ATGCATGCATGCATGCATGC"
#' create_windows(genome, window_size = 4)
#'
#'
#' @export
create_windows <- function(genome, window_size = 1e6) {
  windows <- substring(genome, seq(1, nchar(genome), window_size),
                       seq(window_size, nchar(genome) + window_size, window_size))
  return(windows)
}
