#' Remove Trailing Characters After Space in "Chromosome" Column
#'
#' This function removes everything after the first space in the "Chromosome" column of a data frame.
#' It is useful for standardizing chromosome labels that may have additional text or metadata after a space.
#'
#' @param data A data frame that contains a column named "Chromosome".
#'
#' @return A data frame with modified "Chromosome" values, retaining only the portion of each entry before the first space.
#'
#' @examples
#' # Example usage
#' data <- data.frame(Chromosome = c("Chr1 scaffold", "Chr2 haplotype", "Chr3 variant"))
#' remove_trailing(data)
#'
#' @export
remove_trailing <- function(data) {
  # Function code
  data$Chromosome <- sub(" .*", "", data$Chromosome)
  return(data)
}
