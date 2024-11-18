#' Rename chromosomes by extracting numeric part
#'
#' This function extracts the numeric part of chromosome names in the `Chromosome` column of a dataset. It removes any non-numeric
#' characters and converts the chromosome names into numeric values. This is useful when chromosome names include both letters
#' and numbers, and only the numeric part is required for analysis.
#'
#' @param data A data frame containing a column named `Chromosome` that lists chromosome names. The extraction is based on the numeric part of the chromosome name.
#'
#' @return A data frame with the `Chromosome` column containing only the numeric part of the chromosome names, converted to numeric type.
#'
#' @examples
#' # Example usage
#' data <- data.frame(Chromosome = c("Chr1", "Chr2", "Chr3", "ChrX"),
#'                    Value = c(100, 200, 300, 400))
#' rename_chr(data)
#'
#' @export
rename_chr <- function(data) {

  data$Chromosome <- str_extract(data$Chromosome, "\\d+")  # Extract only numeric
  data$Chromosome <- as.numeric(data$Chromosome)  # make numeric

  return(data)
}
