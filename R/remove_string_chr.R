#' Remove Specific String from "Chromosome" Column
#'
#' This function removes a specified string from the "Chromosome" column of a data frame. 
#' It is helpful for cleaning or standardizing chromosome labels by eliminating unnecessary substrings.
#'
#' @param data A data frame containing a column named "Chromosome".
#' @param remove_string A character string to be removed from the "Chromosome" column. 
#'        Default is "_RagTag".
#'
#' @return A data frame with the specified string removed from the "Chromosome" column.
#'
#' @examples
#' # Example usage
#' data <- data.frame(Chromosome = c("Chr1_RagTag", "Chr2_RagTag", "Chr3"))
#' remove_string_chr(data, remove_string = "_RagTag")
#'
#' @export
remove_string_chr <- function(data, remove_string =  "_RagTag"){
  
  data$Chromosome <- sub(remove_string, "", data$Chromosome)
  
}
