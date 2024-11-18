#' Remove Leading Zeros from Chromosome Numbers
#'
#' This function removes leading zeros from chromosome numbers in the range 01-09 within a specified chromosome prefix.
#'
#' @param data A data frame containing a column named `Chromosome` with chromosome names as strings.
#' @param chr_string A character string representing the prefix for chromosome names (e.g., "Chr"). Default is "Chr".
#'
#' @return A data frame with the `Chromosome` column modified to remove leading zeros from chromosome numbers 01-09.
#'
#' @details
#' The function identifies and removes leading zeros from chromosome identifiers that match the specified prefix in the format `Chr01`, `Chr02`, ..., `Chr09`.
#' This adjustment ensures chromosome numbers are consistently formatted for downstream analysis.
#'
#' @examples
#' # Example usage
#' data <- data.frame(Chromosome = c("Chr01", "Chr02", "Chr10"))
#' remove_lead_0s(data, chr_string = "Chr")
#' # Output will be a data frame with Chromosome values: "Chr1", "Chr2", "Chr10"
#'
#' @export
remove_lead_0s <- function(data, chr_string = "Chr") {

  # Remove leading zero from chromosome numbers 01-09
  data$Chromosome <- sub(paste0("^", chr_string, "0(\\d)$"), paste0(chr_string, "\\1"), data$Chromosome)

  return(data)
}
