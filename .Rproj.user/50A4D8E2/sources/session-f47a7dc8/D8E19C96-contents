#' Select chromosomes from a dataset
#'
#' This function filters a dataset to retain only rows corresponding to chromosomes, based on the presence of a specific string
#' (default is "Chr") in the `Chromosome` column. It removes any extraneous scaffolds or other non-chromosomal data.
#'
#' @param data A data frame containing a column named `Chromosome` that lists the chromosome names.
#' @param chr_string A character string used to filter the chromosome names. Rows with chromosome names containing this string
#'        will be retained. Default is "Chr".
#'
#' @return A data frame containing only the rows where the `Chromosome` column contains the specified `chr_string`.
#'
#' @examples
#' # Example usage
#' data <- data.frame(Chromosome = c("Chr1", "Chr2", "Scaffold1", "Chr3"),
#'                    Value = c(100, 200, 300, 400))
#' select_chr(data, chr_string = "Chr")
#'
#' @export
select_chr <- function(data, chr_string = "Chr"){

  # Remove extraneous scaffolds
  data <- data[grepl(chr_string, data$Chromosome),]

  return(data)
}
