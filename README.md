ggenomics
================
Matthew Davis
2024-11-25

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Functions](#functions)
- [Usage](#usage)
- [Arguments](#arguments)
- [Tutorials](#tutorials)
- [Legacy Functions](#legacy)
- [Getting Help](#help)
- [Contribution](#contribution)
- [License](#license)

## Introduction

`ggenomics` is an R package that provides data visualizations using
ggplot2. It offers functions to dynamically plot genomes for exploratory
data analysis. `ggenomics` aims to utilize ggplot syntax to provide
base-level genomic plots that can later be customized by the user.

`ggenomics` was designed to simplify genomic data visualization using
`ggplot2`. `ggenomics` focuses on:

- Seamless integration with the `ggplot2` ecosystem.
- Dynamic plotting for exploratory analysis.
- Support for large genomic datasets.

## Installation

### Installing Dependencies

To install required dependencies, you can use the following code:

``` r
install.packages(c("data.table", "tidyverse", "scales", "pbapply"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("Biostrings")
```

### Installing `ggenomics`

You can install `ggenomics` from GitHub using the following command:

``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("matthewwdavis/ggenomics")
```

### Verifying Installation

After installing `ggenomics`, load it into R and check the version:

``` r
library(ggenomics)
packageVersion("ggenomics")
```

## Functions

The functions in `ggenomics` create specifically structured data frames
for plotting. `ggenomics` is expected to continually evolve, with more
functions for analysis and plotting added overtime.

Below are some current functions in `ggenomics` and a very brief
description. For more information see [Usage](#usage), and for in-depth
examples from start-to-finish see [Tutorials](#tutorials).

**Current `ggenomics` functions:**  

- `ggread_fasta()` reads in fasta files.
- `telomere_plotting_table()` generates data in a format necessary for
  `geom_telplot()`.
- `ggenom()` initializes a ggplot2 object with `ggenomics` specific
  mapping options.
- `geom_telplot()` creates a plot of chromosomes with telomeric
  sequences marked by size.
- `create_window_fasta()` creates windows from a fasta file read in with
  `ggread_fasta()` or `readDNAStringSet()`.
- `sliding_window_table()` creates sliding windows from a table with
  columns CHROM and POS.

There are many functions that are rarely used on their own and are
instead used to facilitate other, larger functions within `ggenomics`.
Those functions will not receive in-depth documentation, but they are
available as separate functions for the user regardless. The code behind
these functions can be viewed in R with `View(function_name)`.

## Usage

`ggenomics` has two main functionalities: data analysis and plotting.
The data analysis tools are set up to be used with the plotting
functions. A typical workflow will use a specific data analysis tool to
generate a data set with specific formatting. This data set will then be
incorporated into the respective plotting functions.

With the goal of replicating `ggplot2` syntax, `ggenomics` uses a
wrapper function, `ggenom()` to read in data created by other functions.
`ggenom()` has `ggenomics` specific mapping options for plotting with
`ggnomics` geoms. Plotting functions will be added as with `ggplot2`
syntax, attaching geoms to `ggenom()` with a `+`.

This section will describe the functions and cover some basic examples
for each function. The [Tutorials](#tutorials) section has more in-depth
start-to-finish information on usage, and further information about the
arguments available in each function can be found in the
[Arguments](#arguments) section.

`ggread_fasta()`

- This function reads a fasta file into R. It is a wrapper for
  `readDNAStringSet()` and creates a DNAStringSet object. `ggenomics`
  functions that use fasta files will use this object, so using this
  function or `readDNAStringSet()` for fasta files is necessary.

``` r
genome <- ggread_fasta("path/to/fasta")
```

`telomere_plotting_table()`

- This function takes a fasta file and telomere string (Default =
  “CCCTAAA”), then looks through the fasta for the occurrence of that
  telomere string in specified windows (Default = 1 mb). The function
  will look for three of the specified telomere string back-to-back to
  minimize detecting the kmer not associated with telomeres. For
  example, if the string is “CCCTAAA”, the windows of the genome will be
  searched for occurrences of “CCCTAAACCCTAAACCCTAAA”. The table is then
  filtered for a minimum number of string occurrences per window
  (Default = 25).

``` r
tel.table <- telomere_plotting_table(genome)
```

`ggenom()`

- This function initializes a `ggplot2` object in R. It is a wrapper for
  `ggplot()` and creates a `ggplot2` object with `ggenomics` specific
  mapping specified by the `plot` argument. The user should specify
  specific `plot` values to use with corresponding `ggenomics` geoms.
  The proper `plot` argument settings will be in the geom specific
  [Usage](#usage) information and in the [Tutorials](#tutorials) for
  creating specific plots. The options for `plot` and the corresponding
  geoms can be seen in [Arguments](#arguments)

``` r
ggenom_object <- ggenom(tel.table, plot = "telplot")
```

`geom_telplot()`

- This function generates a `ggenomics` telomere plot that can be used
  as a base for genome visualization. It should be added to a `ggenom()`
  created object with the proper mapping information specified by the
  `plot` argument. The proper `plot` argument value to be specified in
  `ggenom()` for `geom_telplot` is `plot = "telplot"`. This can be seen
  here and in [Tutorials](#tutorials). Since the plot is `ggplot2`
  based, the user can customize the result however they like with
  `ggplot2`, a concept which is further explored in
  [Tutorials](#tutorials).

``` r
ggenom(tel.table, plot = "telplot") +
  geom_telplot()

# --or-- #

ggenom_object +
  geom_telplot()
```

`create_window_fasta()`

- This function creates windows (Default = 1mb) from a fasta file read
  with `ggread_fasta()` or `readDNAStringSet()`. It extracts sub-strings
  at regularly defined intervals from the fasta file.

``` r
fasta_windows <- create_windows_fasta(genome)
```

`sliding_window_table()`

- This function creates windows (Default = 10kb) with a slide (Default =
  5kb) from data.frames and data.tables with columns named CHROM and
  POS. It will create 3 new columns (WINDOW_START, WINDOW_END,
  POS_WINDOW) and append them to the current data. WINDOW_START is the
  base pair position of the start of the window, WINDOW_END is the base
  pair position of the end of the window, and POS_WINDOW is the base
  pair position of the midpoint of the window. If the user does not want
  the windows to slide, set the `slide_size` argument equal to the
  `window_size` argument. If a column named SOURCE is present, the
  function will automatically take that into account when creating
  windows.

``` r
table_windows <- sliding_window_table(vcf_table)
```

## Arguments

`ggread_fasta()`

- `path_to_fasta`: The directory path to the fasta file of interest to
  read into R. Creates a DNAStringSet object.

`telomere_plotting_table()`

- `genome`: DNAStringSet object of a fasta file. Can be generated with
  `ggread_fasta` or `readDNAStringSet`.
- `chr_names`: A character string indicating the prefix designating
  chromosome names. This is a crucial argument. Default is “Chr”. If the
  chromsomes begin with a number, use “^\d”.
- `string_remove`: A character string to remove from chromosome names.
  Default is “\_RagTag”. If the user does not want to remove strings
  other than the default, this does not need to be changed.
- `tel_start_seq`: A character string representing the telomere
  sequence. Default is the Arabidopsis telomere repeat, “CCCTAAA”.
- `tel_end_seq`: A character string representing the reverese complement
  of the telomere sequence. Default is the Arabidopsis telomere repeat,
  “TTTAGGG”.
- `size_windows`: A numeric value specifying the size of the window to
  search for telomeric sequence within. Default is 1000000 (1mb).
- `min_tel_count`: A numeric value specifying the minimum telomere
  repeat count per window to include in the final table. Default is 25.
- `sample_name`: A character string to include the sample name in the
  table. Default is NULL.

`ggenom()`

- `data`: Data needed necessary for plotting. This is generated from
  other functions, such as `telomere_plotting_table`.

- `mapping`: The column headers for plotting. If using the `plot`
  argument, this can be ignored. It is suggested to use `plot`.

- `plot`: The necessary mapping information for each `ggenomics` style
  plot. This is different for each `ggenomics` geom specified. This is
  the suggested usage over mapping. Options include:

  - “telplot” to be used with `geom_telplot()`

`geom_telplot()`

- `chr_color`: A character string specifying the color of the plotted
  chromosomes or sequences. Default is “\#F8766D”.
- `chr_size`: A numeric value specifying the width of the plotted
  chromosomes or sequences. Default is 6.
- `tel_color`: A character string specifying the color of the plotted
  telomeres. Default is “black”.
- `tel_shape`: A numeric value specifying the shape of the plotted
  telomeres. Default is 16.
- `legend_title`: A character string specifying the title of the legend.
  Default is “Telomere Size (bp)”.
- `text_size`: A numeric value specifying the base size for plot text
  like axis labels and legends. Default is 6.
- `plot_title`: A character string specifying the title of the plot
  Default is NULL.
- `x_axis_title`: A character string specifying the title of the x-axis
  Default is NULL.
- `y_axis_title`: A character string specifying the title of the y-axis
  Default is “Chromosome Length”.

`create_window_fasta()`

- `genome`: DNAStringSet object of a fasta file. Can be generated with
  `ggread_fasta()` or `readDNAStringSet()`.
- `window_size`: A numeric value specifying the size of the window.
  Default is 1000000 (1mb).

`sliding_window_table()`

- `mut_table`: A data.frame or data.table containing genomic data.
  **Must have columns CHROM and POS**.
- `window_size`: Numeric value specifying size of the window. Default is
  10000 (10kb).
- `slide_size`: Numeric value specifying step size for sliding window.
  If the user does not want slide, the step size should equal the window
  size. Default is 5000 (5kb).

## Tutorials

The following examples are more in-depth than what is found in
[Usage](#usage) and meant to walk the user through using the package
from start-to-finish, data download to plotting. The sub-headers define
different end goals. Publicly available data is used so that the users’
results can be compared here to make sure everything is operating
correctly.

### Creating telomere plots with geom_telplot()

Downloading an example fasta file (Arabidopsis TAIR10):

``` r
download.file("https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz", destfile = "./arabidopsis_tair10.fasta.gz", mode = "wb")
```

The first step is to load the library

``` r
library(ggenomics)
```

After loading the library, read in the fasta file with `ggread_fasta()`

Read in the example fasta to use for `ggenomics`:  
- This creates a DNAStringSet object of a fasta file of interest for
downstream analysis.

``` r
genome <- ggread_fasta("./arabidopsis_tair10.fasta.gz")
```

Next the user should use a data analysis function :  
- In this example, the function creates a table with telomere counts.

``` r
telo.table <- telomere_plotting_table(genome, chr_names = "^\\d")
# "^\\d" is used here to specify that the chromosome names begin with a number, as we are not interested in plotting the plasmid genomes.

print(telo.table)
```

    ##    Chromosome   Length Forward_Counts Reverse_Counts begin_telo_bp end_telo_bp
    ##        <fctr>    <int>          <int>          <int>         <num>       <num>
    ## 1:          1 30427671            270              1          5670          21
    ## 2:          2 19698289              0             49             0        1029
    ## 3:          3 23459830             37              0           777           0
    ## 4:          4 18585056             52              0          1092           0
    ## 5:          5 26975502             NA             NA            NA          NA
    ##    begin_telo_start begin_telo_end end_telo_start end_telo_end total_telo_bp
    ##               <num>          <num>          <num>        <int>         <num>
    ## 1:                0           5670       30427650     30427671          5691
    ## 2:                0              0       19697260     19698289          1029
    ## 3:                0            777       23459830     23459830           777
    ## 4:                0           1092       18585056     18585056          1092
    ## 5:                0             NA             NA     26975502            NA
    ##    normalized_total_telo_size
    ##                         <num>
    ## 1:               4.776479e-05
    ## 2:               8.636438e-06
    ## 3:               6.521392e-06
    ## 4:               9.165199e-06
    ## 5:                         NA

**NOTE:** It is always a good idea to inspect the table and ensure you
are seeing what is expected. In this case, Chromosome 5 had no detected
telomeric repeat, and so it has NA values.

Then, the user can utilize the `ggenom()` function paired with the
`geom_telplot()` function to create a telomere plot:  
- If the user wants to create a telomere plot, the set
`plot = "telplot"`. The options for possible `plot` argument values can
be found in [Arguments](#arguments).

``` r
ggenom(telo.table, plot = "telplot") +
  geom_telplot()
```

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

There are some arguments within `geom_telplot()` to specify shape and
color:

``` r
ggenom(telo.table, plot = "telplot") +
  geom_telplot(chr_color = "bisque2", tel_color = "darkgreen", tel_shape = 18)
```

![](README_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Since all plots are `ggplot2` based, they can be edited and adjusted
like any ggplot:  
- The user can add adjustments with `+`, just like in `ggplot2`.

``` r
ggenom(telo.table, plot = "telplot") +
  geom_telplot(chr_color = "bisque2", tel_color = "darkgreen", tel_shape = 18) +
  scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "Mb")) +
  labs(y = "Sequence Length", x = "Chromosome", size = "Telomere Size", title = "ggenomics Telomere Plot") +
  theme_classic(base_size = 6) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.2, "cm"),
        plot.title = element_text(hjust = 0.5, face = "bold"))
```

![](README_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

## Legacy Functions

**Disclaimer:** These functions may not be actively maintained, and
users should use updated alternatives when possible.

These functions were originally from the `ggideo` package. While that
package has been archived, the functions will continue to exist in
`ggenomics`, albeit with little continuous upkeep.

Below is an example of how to use these legacy functions:

- `ggideo()` is used to plot telomere plots of primary assemblies.

``` r
library(ggenomics)

# Generate data and plot, stored as a list
genome.plot <- ggideo("./arabidopsis_tair10.fasta.gz", chr_names = "^\\d")
```

- Print the table.

``` r
genome.plot$genomic.table
```

- Print the plot.

``` r
genome.plot$ideogram
```

![](README_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

- `ggideo_diploid()` is used to plot telomere plots of haplotype phased
  diploid assemblies. The haplotypes can be two separate fasta files, or
  a fasta file with both haplotypes present. The haplotypes should be
  identified with “\_hap1” and “\_hap2”.
- Example of the two separate fasta files.

``` r
library(ggenomics)

# Generate data and plot, stored as a list
genome.plot <- ggideo_diploid("./genome_hap1.fasta.gz", "./genome_hap2.fasta.gz")
```

    ## Joining with `by = join_by(Chromosome, Length, Forward_Counts, Reverse_Counts,
    ## begin_telo_bp, end_telo_bp, begin_telo_start, begin_telo_end, end_telo_start,
    ## end_telo_end, total_telo_bp, normalized_total_telo_size, Hap)`

- Print the table.

``` r
genome.plot$genomic.table
```

- Print the plot.

``` r
genome.plot$ideogram
```

![](README_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

- Example of usage with both haplotypes in one combined fasta file.

``` r
library(ggenomics)

# Generate data and plot, stored as a list
genome.plot <- ggideo_diploid(combined_hap_fasta = "./genome_combohap.fasta.gz",
                              string_remove = "_hap\\d_RagTag")
```

    ## Joining with `by = join_by(Chromosome, Length, Hap, Forward_Counts,
    ## Reverse_Counts, begin_telo_bp, end_telo_bp, begin_telo_start, begin_telo_end,
    ## end_telo_start, end_telo_end, total_telo_bp, normalized_total_telo_size)`

- Print the table.

``` r
genome.plot$genomic.table
```

- Print the plot.

``` r
genome.plot$ideogram
```

![](README_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

## Getting Help

If you encounter any issues or have questions, you can:

- Check the detailed package documentation using `?function_name` in R
- Submit an issue on the [GitHub
  repository](https://github.com/matthewwdavis/ggenomics/issues)

## Contribution

Contributions are welcome! To contribute:

- Fork the repository on GitHub.
- Make your changes in a new branch.
- Submit a pull request with a detailed description of your changes.

## License

This package is licensed under the MIT License. See the
[LICENSE](https://github.com/matthewwdavis/ggenomics/blob/main/LICENSE)
file for details.
