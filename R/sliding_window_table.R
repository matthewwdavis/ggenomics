#' Generate Sliding Windows From Table
#'
#' This function calculates sliding windows for a mutation table, enabling downstream analysis of genomic regions.
#' It computes window start and end positions as well as the midpoint for each window, based on `CHROM` and `POS` columns.
#'
#' **Important:** The input data frame must contain columns named `CHROM` (chromosome) and `POS` (position).
#'
#' @param mut_table A data frame containing mutation data. **Must include columns named `CHROM` and `POS`.**
#' @param window_size A numeric value specifying the size of the sliding window (in base pairs). Default is `10000` (10 kb).
#' @param slide_size A numeric value specifying the step size for sliding the window (in base pairs). Default is `5000` (5 kb).
#'
#' @return A data frame with the original data and additional columns:
#' \describe{
#'   \item{WINDOW_START}{Start position of the sliding window.}
#'   \item{WINDOW_END}{End position of the sliding window.}
#'   \item{POS_WINDOW}{Midpoint of the sliding window.}
#' }
#'
#' @examples
#' # Example input data
#' mut_table <- data.frame(
#'   CHROM = c("Chr1", "Chr1", "Chr1", "Chr2", "Chr2"),
#'   POS = c(100, 5000, 15000, 200, 8000)
#' )
#'
#' # Generate sliding window table
#' sliding_window_table(mut_table, window_size = 10000, slide_size = 5000)
#'
#' @importFrom dplyr arrange group_by mutate ungroup
#' @export
sliding_window_table <- function(mut_table, window_size = 10000, slide_size = 5000) {
  # Check if the SOURCE column exists
  if ("SOURCE" %in% colnames(mut_table)) {
    # Include SOURCE in arrange and group_by
    table_w_windows <- mut_table %>%
      arrange(SOURCE,CHROM, POS) %>%
      group_by(SOURCE, CHROM) %>%
      mutate(
        WINDOW_START = ((POS - 1) %/% slide_size) * slide_size + 1,
        WINDOW_END = WINDOW_START + window_size - 1,
        POS_WINDOW = ((WINDOW_START + WINDOW_END) / 2)
      ) %>%
      ungroup()
  } else {
    # Original behavior
    table_w_windows <- mut_table %>%
      arrange(CHROM, POS) %>%
      group_by(CHROM) %>%
      mutate(
        WINDOW_START = ((POS - 1) %/% slide_size) * slide_size + 1,
        WINDOW_END = WINDOW_START + window_size - 1,
        POS_WINDOW = ((WINDOW_START + WINDOW_END) / 2)
      ) %>%
      ungroup()
  }
  
  return(table_w_windows)
}

