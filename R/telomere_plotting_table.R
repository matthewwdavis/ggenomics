telomere_plotting_table <- function(genome, chr_names = "Chr", string_remove = "_RagTag", tel_start_seq = "CCCTAAA", tel_end_seq = "TTTAGGG",
                          size_windows = 1e6, min_tel_count = 25, sample_name = NULL){
  # Create table of contigs, chromosomes, and lengths
  length.table <- data.table(Chromosome = names(genome), Length = width(genome))
  
  # Filter for only chromosomes, based on starting string
  length.table <- select_chr(length.table, chr_string = chr_names)
  
  # Extract the size of the genome
  genome.size <- sum(length.table$Length)
  
  # Count telomeric sequence repeat
  tel_count.table <- telomere_repeat_number(fasta = genome, window = size_windows, tel_start = tel_start_seq, tel_end = tel_end_seq)
  
  # Filter to maintain telomeric counts over a certain threshold
  tel.table <- tel_count.table %>%
    filter(Forward_Counts >= min_tel_count | Reverse_Counts >= min_tel_count)
  
  # Create the larger table necessary for plotting
  plotting.table <- genome_table(length.table, tel.table, name = sample_name, genome_size = genome.size)
  
  # Remove string from names
  plotting.table$Chromosome <- remove_string_chr(plotting.table, remove_string = string_remove)
  
  # Remove leading 0s for proper ordering and plotting
  plotting.table <- remove_lead_0s(plotting.table, chr_string = chr_names)
  
  # Remove trailing strings. These is defined by a space
  plotting.table <- remove_trailing(plotting.table)
  
  # Set levels so that chromosomes are plotted in the proper order by number
  plotting.table$Chromosome <- factor(plotting.table$Chromosome,
                                      levels = unique(plotting.table$Chromosome)[
                                        order(sapply(unique(plotting.table$Chromosome), function(x) {
                                          # If the chromosome value is purely numeric, keep it as is
                                          if (grepl("^\\d+$", x)) {
                                            return(as.numeric(x))  # Convert numeric strings to numbers
                                          } else {
                                            # For non-numeric chromosomes, remove everything after the space
                                            return(as.character(gsub(" .*", "", x)))
                                          }
                                        }))
                                      ])
  
  return(plotting.table)
}
