# HTGQC

## HTGQC is an R package for quality check of HTG EdgeSeq datasets.

To install, please install the "devtools" R package, then run

`devtools::install_github("LodovicoTerzi/HTGQC")`

The package has two main functions:

`readHTG` reads the unformatted Excel file as it outputs from the HTG Edge-seq machines.
  
  The function prepares the unformatted file for the quality control.
  
`qualityCheck` performs the quality check of the data
  
  The function takes as input (i) the output from readHTG function, or (ii) a data frame containing gene names as row names and sample IDs as column names.
  The user can specify the number of positive and negative controls, their order, and the path where to write the output files.
  The function writes a folder with one plot and one table for the positive and negative control and, if specified, a cleaned table without controls, and a normalised table (log2 cpm).

An R shiny application has been developed for researchers with no knowledge of R. Please visit `https://lodovico.shinyapps.io/shinyHTGQC` for details.
