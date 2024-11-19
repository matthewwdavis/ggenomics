#' Extract Haplotype Information from "Chromosome" Column
#'
#' This function extracts haplotype identifiers from the "Chromosome" column of a data frame using a specified pattern 
#' and creates a new column named "Hap" to store the extracted haplotype information.
#'
#' @param data A data frame containing a column named "Chromosome".
#' @param hap_string A character string representing the regular expression pattern used to identify haplotypes in the "Chromosome" column. 
#'        Default is `"hap\\d+"` (matches substrings like "hap1", "hap2", etc.).
#'
#' @return A data frame with an additional column named "Hap" containing the extracted haplotype information.
#'         If no match is found, the "Hap" column will contain `NA` for those rows.
#'
#' @examples
#' # Example usage
#' data <- data.frame(Chromosome = c("hap1_Chr1", "hap2_Chr2", "hap3_Chr3"))
#' parse_hap_combo(data, hap_string = "hap\\d+")
#'
#' @importFrom dplyr mutate
#' @importFrom stringr str_extract
#' @export
parse_hap_combo <- function(data, hap_string = "hap\\d+") {
  # Extract haplotype information and add as a new column
  data <- data %>%
    
    mutate(Hap = str_extract(Chromosome, hap_string))
  
  return(data)
}
