#' Count occurrences of a sequence within a window
#'
#' This function counts the number of occurrences of a given sequence within a specified window.
#' It uses the `countPattern` function to count the exact occurrences of the sequence in the window.
#'
#' @param window A character string representing the window of DNA or RNA sequence to search within.
#' @param sequence A character string representing the sequence to count within the window.
#'
#' @return An integer value representing the total count of occurrences of the sequence within the window.
#'
#' @examples
#' # Example usage
#' window <- "AGCTAGCTAGCT"
#' sequence <- "AGC"
#' count_sequence(window, sequence)
#'
#' @seealso \code{\link[Biostrings]{countPattern}} for more details on pattern matching.
#'
#' @export
count_sequence <- function(window, sequence) {
  count <- sum(countPattern(sequence, window))
  return(count)
}
