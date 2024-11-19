#' Generate a diploid haplotype phased ideogram plot from two fasta file
#'
#' This function creates an ideogram plot based on a genome sequence provided in a FASTA file. It reads the genome,
#' filters for specified chromosomes, calculates telomeric repeat counts, and visualizes the chromosome lengths with telomere regions marked.
#' The output includes a data table and a `ggplot` generated ideogram plot.
#'
#' @param path_hap1_fasta A character string specifying the file path to the fasta file containing the haplotype 1 genome sequence.
#' @param path_hap2_fasta A character string specifying the file path to the fasta file containing the haplotype 2 genome sequence.
#' @param chr_names A character string specifying the pattern to identify chromosome names in the data. Default is "Chr".
#' @param tel_start_seq A character string representing the telomeric repeat sequence at the start of the chromosome. Default is "CCCTAAA".
#' @param tel_end_seq A character string representing the telomeric repeat sequence at the end of the chromosome. Default is "TTTAGGG".
#' @param size_windows A numeric value specifying the window size (in base pairs) for counting telomeric repeats. Default is `1e6` (1 Mb).
#' @param min_tel_count A numeric value representing the minimum number of telomeric repeats required to retain telomeric counts in the plot. Default is 25.
#' @param sample_name A character string specifying the name of the sample for labeling purposes. Default is `NULL`.
#' @param title_plot A character string specifying the title of the plot. Default is `NULL`.
#' @param title_x_axis A character string specifying the title of the x axis. Default is `NULL`.
#' @param title_y_axis A character string specifying the title of the y axis. Default is "Chromosome Length".
#' @param title_legend A character string specifying the title of the legend. Default is "Telomere Presence".
#' @param color_hap1 A character string specifying the color of the chromosome segments. Default is "dodgerblue2".
#' @param color_hap2 A character string specifying the color of the chromosome segments. Default is "orangered".
#' @param name_hap1 A character string specifying the name of the haplotype. Default is "Haplotype 1".
#' @param name_hap2 A character string specifying the name of the haplotype. Default is "Haplotype 2".
#' @param size_chr A numeric value specifying the size (linewidth) of the chromosome segments. Default is 6.
#' @param color_tel A character string specifying the color of the telomere points. Default is "black".
#' @param shape_tel A numeric value specifying the shape of the telomere points. Default is 16 (filled circle).
#' @param scale_y A numeric value for scaling the y-axis. Default is `1e-6` to convert lengths to megabases (Mb).
#' @param suffix_y_scale A character string to append to the y-axis labels, usually a unit like "Mb". Default is "Mb".
#' @param pos_legend A character string specifying the position of the legend. Default is "bottom".
#' @param size_legend A numeric value specifying the size of the legend keys. Default is 0.25.
#' @param size_text A numeric value specifying the size of text for the entire figure. Adjusts `base_size` ggplot function. Default is 6.
#' @return A list containing:
#' \describe{
#'   \item{genomic.table}{A data frame with chromosome, length, and telomere information.}
#'   \item{ideogram}{A `ggplot` object representing the ideogram plot.}
#' }
#'
#' @examples
#' # Example usage
#' # ggideo_diploid(path_hap1_fasta = "path/to/genome.fasta", path_hap2_fasta = "path/to/genome.fasta", sample_name = "Sample_1", title_plot = "Genome Ideogram")
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom data.table data.table
#' @importFrom ggplot2 ggplot aes geom_segment geom_point labs theme element_text scale_y_continuous
#' @importFrom scales label_number
#' @importFrom dplyr filter mutate left_join select
#'
#' @export
# ggideo_diploid <- function(path_hap1_fasta, path_hap2_fasta, chr_names = "Chr", tel_start_seq = "CCCTAAA", tel_end_seq = "TTTAGGG",
#                    size_windows = 1e6, min_tel_count = 25, sample_name = NULL, title_plot = NULL,
#                    title_x_axis = NULL, title_y_axis = "Chromosome Length", title_legend = "Telomere Size (bp)",
#                    color_hap1 = "dodgerblue2", color_hap2 = "orangered", name_hap1 = "Haplotype 1",
#                    name_hap2 = "Haplotype 2", size_chr = 6, color_tel = "black", shape_tel = 16,
#                    scale_y = 1e-6, suffix_y_scale = "Mb", pos_legend = "bottom", size_legend = 0.25, size_text = 6){
#   
#   # Read in fasta files
#   hap1.genome <- readDNAStringSet(path_hap1_fasta)
#   hap2.genome <- readDNAStringSet(path_hap2_fasta)
#   
#   # Create table of contigs, chromosomes, and lengths
#   hap1_length.table <- data.table(Chromosome = names(hap1.genome), Length = width(hap1.genome))
#   hap2_length.table <- data.table(Chromosome = names(hap2.genome), Length = width(hap2.genome))
#   
#   # Filter for only chromosomes, based on starting string
#   hap1_length.table <- select_chr(hap1_length.table, chr_string = chr_names)
#   hap2_length.table <- select_chr(hap2_length.table, chr_string = chr_names)
#   
#   # Extract the size of the genome
#   hap1_genome.size <- sum(hap1_length.table$Length)
#   hap2_genome.size <- sum(hap2_length.table$Length)
#   
#   # Count telomeric sequence repeat
#   hap1_tel_count.table <- telomere_repeat_number(fasta = hap1.genome, window = size_windows, tel_start = tel_start_seq, tel_end = tel_end_seq)
#   hap2_tel_count.table <- telomere_repeat_number(fasta = hap2.genome, window = size_windows, tel_start = tel_start_seq, tel_end = tel_end_seq)
#   
#   # Filter to maintain telomeric counts over a certain threshold
#   hap1_tel.table <- hap1_tel_count.table %>%
#     filter(CCCTAAA_Counts >= min_tel_count | TTTAGGG_Counts >= min_tel_count)
#   
#   hap2_tel.table <- hap2_tel_count.table %>%
#     filter(CCCTAAA_Counts >= min_tel_count | TTTAGGG_Counts >= min_tel_count)
#   
#   # Create the larger table necessary for plotting
#   hap1_plotting.table <- genome_table(hap1_length.table, hap1_tel.table, name = sample_name, genome_size = hap1_genome.size)
#   hap2_plotting.table <- genome_table(hap2_length.table, hap2_tel.table, name = sample_name, genome_size = hap2_genome.size)
#   
#   # Assign haplotypes to individuals
#   hap1_plotting.table$Hap <- "hap1"
#   hap2_plotting.table$Hap <- "hap2"
#   diploid_plotting.table <- full_join(hap1_plotting.table, hap2_plotting.table)
#   
#   # Remove leading 0s for proper ordering and plotting
#   diploid_plotting.table <- remove_lead_0s(diploid_plotting.table, chr_string = chr_names)
#   
#   # Remove trailing strings. These is defined by a space
#   diploid_plotting.table <- remove_trailing(diploid_plotting.table)
#   
#   # Set levels so that chromosomes are plotted in the proper order by number
#   diploid_plotting.table$Chromosome <- factor(diploid_plotting.table$Chromosome,
#                                       levels = unique(diploid_plotting.table$Chromosome)[order(as.numeric(gsub(chr_names, "", unique(diploid_plotting.table$Chromosome))))])
#   
#   # plot the ideogram
#   graphic <- ideogram_diploid(genome.table = diploid_plotting.table, plot_title = title_plot, x_axis_title = title_x_axis,
#                       y_axis_title = title_y_axis, legend_title = title_legend, hap1_color = color_hap1,
#                       hap2_color = color_hap2, hap1_name = name_hap1, hap2_name = name_hap2,
#                       chr_size = size_chr, tel_color = color_tel, tel_shape = shape_tel,
#                       y_scale = scale_y, y_scale_suffix = suffix_y_scale, legend_pos = pos_legend,
#                       legend_size = size_legend, text_size = size_text)
#   
#   return(list(genomic.table = diploid_plotting.table, ideogram = graphic))
# }

ggideo_diploid <- function(path_hap1_fasta, path_hap2_fasta, combined_hap_fasta = NULL, chr_names = "Chr", 
                           string_remove = "_RagTag", string_hap = "hap\\d+",tel_start_seq = "CCCTAAA", tel_end_seq = "TTTAGGG",
                           size_windows = 1e6, min_tel_count = 25, sample_name = NULL, title_plot = NULL,
                           title_x_axis = NULL, title_y_axis = "Chromosome Length", title_legend = "Telomere Size (bp)",
                           color_hap1 = "dodgerblue2", color_hap2 = "orangered", name_hap1 = "Haplotype 1",
                           name_hap2 = "Haplotype 2", size_chr = 6, color_tel = "black", shape_tel = 16,
                           scale_y = 1e-6, suffix_y_scale = "Mb", pos_legend = "bottom", size_legend = 0.25, size_text = 6){
  
  if (!is.null(combined_hap_fasta)) {
    
    # Read in combined fasta file
    combined.genome <- readDNAStringSet(combined_hap_fasta)
    
    # Create table of contigs, chromosomes, and lengths
    combined_length.table <- data.table(Chromosome = names(combined.genome), Length = width(combined.genome))
    
    # Filter for only chromosomes, based on starting string
    combined_length.table <- select_chr(combined_length.table, chr_string = chr_names)
    
    # Remove string from names
    combined_length.table <- remove_string_chr(combined_length.table, remove_string = string_remove)
    
    # Extract the size of the genome
    combined_genome.size <- sum(combined_length.table$Length)
    
    # Count telomeric sequence repeat
    combined_tel_count.table <- telomere_repeat_number(fasta = combined.genome, window = size_windows, tel_start = tel_start_seq, tel_end = tel_end_seq)
    
    # Filter to maintain telomeric counts over a certain threshold
    combined_tel.table <- combined_tel_count.table %>%
      filter(CCCTAAA_Counts >= min_tel_count | TTTAGGG_Counts >= min_tel_count)
    
    # Create the larger table necessary for plotting
    combined_plotting.table <- genome_table(combined_length.table, combined_tel.table, name = sample_name, genome_size = combined_genome.size)
    
    # Assign haplotypes for plotting table
    combined_plotting.table <- parse_hap_combo(combined_plotting.table, hap_string = string_hap)
    
    # Remove leading 0s for proper ordering and plotting
    combined_plotting.table <- remove_lead_0s(combined_plotting.table, chr_string = chr_names)
    
    # Remove trailing strings. These are defined by a space
    combined_plotting.table <- remove_trailing(combined_plotting.table)
    
    # Set levels so that chromosomes are plotted in the proper order by number
    combined_plotting.table$Chromosome <- factor(combined_plotting.table$Chromosome,
                                                 levels = unique(combined_plotting.table$Chromosome)[order(as.numeric(gsub(chr_names, "", unique(combined_plotting.table$Chromosome))))])
    
    # Plot the ideogram
    graphic <- ideogram_diploid(genome.table = combined_plotting.table, plot_title = title_plot, x_axis_title = title_x_axis,
                                y_axis_title = title_y_axis, legend_title = title_legend, hap1_color = color_hap1,
                                hap2_color = color_hap2, hap1_name = name_hap1, hap2_name = name_hap2,
                                chr_size = size_chr, tel_color = color_tel, tel_shape = shape_tel,
                                y_scale = scale_y, y_scale_suffix = suffix_y_scale, legend_pos = pos_legend,
                                legend_size = size_legend, text_size = size_text)
    
    return(list(genomic.table = combined_plotting.table, ideogram = graphic))
    
  } else {
    
    # Read in fasta files
    hap1.genome <- readDNAStringSet(path_hap1_fasta)
    hap2.genome <- readDNAStringSet(path_hap2_fasta)
    
    # Create table of contigs, chromosomes, and lengths
    hap1_length.table <- data.table(Chromosome = names(hap1.genome), Length = width(hap1.genome))
    hap2_length.table <- data.table(Chromosome = names(hap2.genome), Length = width(hap2.genome))
    
    # Filter for only chromosomes, based on starting string
    hap1_length.table <- select_chr(hap1_length.table, chr_string = chr_names)
    hap2_length.table <- select_chr(hap2_length.table, chr_string = chr_names)
    
    # Remove string from names
    hap1_length.table <- remove_string_chr(hap1_length.table, remove_string = string_remove)
    hap2_length.table <- remove_string_chr(hap2_length.table, remove_string = string_remove)
    
    # Extract the size of the genome
    hap1_genome.size <- sum(hap1_length.table$Length)
    hap2_genome.size <- sum(hap2_length.table$Length)
    
    # Count telomeric sequence repeat
    hap1_tel_count.table <- telomere_repeat_number(fasta = hap1.genome, window = size_windows, tel_start = tel_start_seq, tel_end = tel_end_seq)
    hap2_tel_count.table <- telomere_repeat_number(fasta = hap2.genome, window = size_windows, tel_start = tel_start_seq, tel_end = tel_end_seq)
    
    # Filter to maintain telomeric counts over a certain threshold
    hap1_tel.table <- hap1_tel_count.table %>%
      filter(CCCTAAA_Counts >= min_tel_count | TTTAGGG_Counts >= min_tel_count)
    
    hap2_tel.table <- hap2_tel_count.table %>%
      filter(CCCTAAA_Counts >= min_tel_count | TTTAGGG_Counts >= min_tel_count)
    
    # Create the larger table necessary for plotting
    hap1_plotting.table <- genome_table(hap1_length.table, hap1_tel.table, name = sample_name, genome_size = hap1_genome.size)
    hap2_plotting.table <- genome_table(hap2_length.table, hap2_tel.table, name = sample_name, genome_size = hap2_genome.size)
    
    # Assign haplotypes to individuals
    hap1_plotting.table$Hap <- "hap1"
    hap2_plotting.table$Hap <- "hap2"
    diploid_plotting.table <- full_join(hap1_plotting.table, hap2_plotting.table)
    
    # Remove leading 0s for proper ordering and plotting
    diploid_plotting.table <- remove_lead_0s(diploid_plotting.table, chr_string = chr_names)
    
    # Remove trailing strings. These are defined by a space
    diploid_plotting.table <- remove_trailing(diploid_plotting.table)
    
    # Set levels so that chromosomes are plotted in the proper order by number
    diploid_plotting.table$Chromosome <- factor(diploid_plotting.table$Chromosome,
                                                levels = unique(diploid_plotting.table$Chromosome)[order(as.numeric(gsub(chr_names, "", unique(diploid_plotting.table$Chromosome))))])
    
    # Plot the ideogram
    graphic <- ideogram_diploid(genome.table = diploid_plotting.table, plot_title = title_plot, x_axis_title = title_x_axis,
                                y_axis_title = title_y_axis, legend_title = title_legend, hap1_color = color_hap1,
                                hap2_color = color_hap2, hap1_name = name_hap1, hap2_name = name_hap2,
                                chr_size = size_chr, tel_color = color_tel, tel_shape = shape_tel,
                                y_scale = scale_y, y_scale_suffix = suffix_y_scale, legend_pos = pos_legend,
                                legend_size = size_legend, text_size = size_text)
    
    return(list(genomic.table = diploid_plotting.table, ideogram = graphic))
  }
}
