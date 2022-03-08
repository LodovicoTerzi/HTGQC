#'Extract count data from HTG EdgeSeq software's output file
#'
#'Takes in the Excel file output from the HTG EdgeSeq analysis software and performs a quality check, normalisation and annotation.
#'
#'@param data.in
#'The excel file as output from HTG EdgeSeq analysis platform
#'@param path.out
#'The output path where a new folder containing the results will be saved (optional).
#'If not specified, the current working directory is used.
#'
#'@details
#'The quality check is performed on 4 positive control and 4 negative control genes.
#'For the positive controls, samples are considered failures if the percentage of reads allocated to the positive controls is greater than 40%.
#'For the negative controls, samples are considered failures if the deviation from the expected values is greater than twice the standard deviation accross all samples deviations from the expected value.
#'
#'@return
#'A folder contaiing two plots and two summary table, one the for positive and the negative control genes.
#'
#'@export
#'

readHTG <- function(data.in){
  options(stringsAsFactors = F)
  data <- data.in
  data <- data[-c(1:8),]
  row.names(data) <- data[,1]
  colnames(data) <- data[1,]
  data <- data[,-1]
  lib_size <- as.numeric(unlist(data[2,], use.names=F))
  sample_names <- unlist(data[1,], use.names=F)
  data <- data[-c(1,2),]  #what is this?
}
