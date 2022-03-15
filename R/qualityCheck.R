#'Plots and summary of Quality Check for an HTG EdgeSeq run
#'
#'Takes in the Excel file output from the HTG EdgeSeq analysis software and performs a quality check, normalisation and annotation.
#'
#'
#'@param data.in
#'A dataframe of counts. Column names should correspond to sample names, and row names should contain gene names, including positive and negative controls
#'
#'@param TableOut
#'Whether intermediate files will be written.
#'If TRUE, the script will write the clean dataset contaning only the genes of interest, and normalised data.
#'
#'@param PosContrNo
#'Number of positive controls
#'
#'#'@param NegContrNo
#'Number of negative controls
#'
#'#'@param OrderContr
#'Whether the dataframe contains first the positive or negative control genes, in order
#'
#'
#'@details
#'
#'The quality check is performed on 4 positive control and 4 negative control genes.
#'For the positive controls, samples are considered failures if the percentage of reads allocated to the positive controls is greater than 40%.
#'For the negative controls, samples are considered failures if the deviation from the expected values is greater than twice the standard deviation accross all samples deviations from the expected value.
#'
#'@return
#'A folder containing two plots and two summary table, one the for positive and the negative control genes.
#'If TableOut is TRUE, two table are created with clean data (raw gene expression without control genes), and normalised counts (log2 cpm)
#'
#'@export
#'


qualityCheck <- function(data.in, TableOut=FALSE, PosContrNo=4, NegContrNo=4, OrderContr="negative", path.out="current"){

  # set working directory where file is located
  if (path.out == "current"){
    setwd(getwd())
  }
  else {setwd(path.out)}

  # create output directory and move to that directory
  dir.create("QC")
  setwd("QC")

  data <- data.in
  lib_size <- colSums(data)
  sample_names <- colnames(data)

  ## write tables
  if (TableOut==TRUE){
    data_clean <- data[-c(1:(PosContrNo+PosContrNo)),]
    write.table(data_clean, "HTG_CleanData.txt", sep="\t", quote=F, row.names = T)
    write.table(cpm(data_clean, lib.size=colSums(data_clean), log=T), "HTG_CpmData.txt", sep="\t", quote=F, row.names = T)
  }

  #check argument
  assertthat::assert_that(OrderContr %in% c("positive", "negative"), msg = "OrderContr argument should be either 'positive' or 'negative'")

  ## negative controls ##
  if (OrderContr == "negative") {neg <- matrix(as.numeric(as.matrix(data[1:NegContrNo,])), nrow=4)}
  else {neg <- matrix(as.numeric(as.matrix(data[PosContrNo+1:PosContrNo+NegContrNo,])), nrow=4)}

  #plot
  cpm <- cpm(neg, lib.size=lib_size, log=T)
  mean_cpm <- colMeans(cpm)
  mean_all <- mean(mean_cpm)
  delta_cpm <- mean_cpm - mean_all
  stdev_cpm <- sd(mean_cpm)
  png("QC_neg_plot.png", width = 1000, height = 700)
  y_up <- 3*stdev_cpm+1
  if (sum(delta_cpm > 3*stdev_cpm+1)>0) {y_up <- max(delta_cpm)}
  y_down <- -3*stdev_cpm-1
  if (sum(delta_cpm < -3*stdev_cpm-1)>0) {y_down <- min(delta_cpm)}
  plot(c(1:length(delta_cpm)), delta_cpm, ylim=c(y_down,y_up), col="blue3", pch=19, las=1, xaxt='n',
       ylab="Difference from expected value", xlab=NA)
  abline(v=c(1:length(delta_cpm)), col="grey90")
  points(c(1:length(delta_cpm)), delta_cpm, col="blue3", pch=19, las=1, xaxt='n', ylab="Difference from expected value", xlab=NA)
  axis(1, at=c(1:length(delta_cpm)), labels = sample_names, las=2, cex.axis=0.7)
  abline(2*stdev_cpm,0, lty=2, col="green", lwd=2)
  abline(-2*stdev_cpm,0, lty=2, col="green", lwd=2)
  abline(3*stdev_cpm,0, lty=2, col="red", lwd=2)
  abline(-3*stdev_cpm,0, lty=2, col="red", lwd=2)
  dev.off()
  #percentage of reads allocated to negative control
  perc <- colSums(neg)/lib_size
  res <- rep("OK", ncol(neg))
  res[perc > 0.1] <- "ALERT"
  table <- cbind(sample_names, round(perc*100, digits=3), res)
  colnames(table) <- c("SAMPLE", "PercentageOverLibrary", "QC_result")
  print(cbind(sample_names, round(perc*100, digits=3), res))
  write.table(table, "QC_neg_table.txt", quote=F, row.names=F, col.names=T)


  ## positive controls ##
  if (OrderContr == "positive") {pos <- matrix(as.numeric(as.matrix(data[1:PosContrNo,])), nrow=4)}
  else {pos <- matrix(as.numeric(as.matrix(data[NegContrNo+1:NegContrNo+PosContrNo,])), nrow=4)}

  perc2 <- colSums(pos)/lib_size
  res2 <- rep("OK", ncol(pos))
  res2[perc2 > 0.1] <- "ALERT"
  res2[perc2 > 0.4] <- "FAIL"
  table2 <- cbind(sample_names, round(perc2*100, digits=3), res2)
  colnames(table2) <- c("SAMPLE", "PercentageOverLibrary", "QC_result")
  print(cbind(sample_names, round(perc2*100, digits=3), res2))
  write.table(table2, "QC_pos_table.txt", quote=F, row.names=F)
  png("QC_pos_plot.png", width = 700, height = 700)
  plot(c(1:length(delta_cpm)), perc2, ylim=c(0,1), col="blue3", pch=19, las=1, xaxt='n', yaxt='n',
       ylab="% allocated to positive controls", xlab=NA)
  abline(v=c(1:length(delta_cpm)), col="grey90")
  points(c(1:length(delta_cpm)), perc2, ylim=c(0,1), col="blue3", pch=19, las=1, xaxt='n', yaxt='n',
         ylab="% allocated to positive controls", xlab=NA)
  axis(1, at=c(1:length(delta_cpm)), labels = sample_names, las=2, cex.axis=1)
  axis(2, at=c(0,0.2,0.4,0.6,0.8,1), lab=c(0,20,40,60,80,100), las=2)
  abline(0.4,0, lty=2, col="red", lwd=2)
  dev.off()
}
