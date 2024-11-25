#' Read a FASTA File into a DNAStringSet
#'
#' This function reads a DNA sequence from a FASTA file and returns it as a `DNAStringSet` object.
#'
#' @param path_to_fasta A character string specifying the path to the FASTA file.
#'
#' @return A `DNAStringSet` object containing the DNA sequences from the FASTA file.
#'
#' @details This function serves as a wrapper around the `readDNAStringSet` function from the Biostrings package, 
#' to maintain syntax when loading DNA sequences from FASTA files.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' dna_sequences <- ggread.fasta("path/to/your/fasta_file.fasta")
#' }
#'
#' @seealso [Biostrings::readDNAStringSet()]
#'
#' @import Biostrings
#' @export
ggread_fasta <- function(path_to_fasta) {
  readDNAStringSet(path_to_fasta)
}
