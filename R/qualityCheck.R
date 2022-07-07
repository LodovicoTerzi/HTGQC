#'Plots and summary of Quality Check for an HTG EdgeSeq run
#'
#'
#'Performs a quality check, normalisation and annotation of the HTG EdgeSeq data.
#'
#'
#'@param data.in
#'A dataframe of counts. Column names should correspond to sample names, and row names should contain gene names, including positive and negative controls
#'
#'@param TableOut
#'Whether intermediate files will be written.
#'If TRUE, the script will write the clean dataset contaning only the genes of interest (without controls), and normalised data (log counts per million).
#'
#'@param PosContrNo
#'Number of positive controls
#'
#'#'@param NegContrNo
#'Number of negative controls
#'
#'#'@param OrderContr
#'Whether the dataframe contains first the positive or negative control genes, in order from top to bottom
#'
#'@param path.out
#'Directory where the results will be written. A new folder called QC will be created in the specified directory. If not specified, the results will be written in the current directory.
#'
#'@details
#'
#'The quality check is performed on 4 positive control and 4 negative control genes. The control genes should be the first rows in the dataset, followed by the genes profiled.
#'For the positive controls, samples are considered failures if the percentage of reads allocated to the positive controls is greater than 40%.
#'For the negative controls, samples are considered failures if the deviation from the expected values is greater than twice the standard deviation (see manuscript for details).
#'
#'@return
#'A folder containing two plots and two summary table, one the for positive and the negative control genes.
#'If TableOut is TRUE, two additional tables are created with clean data (raw gene expression without control genes), and normalised counts (log2 cpm)
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
  if (OrderContr == "negative") {neg <- matrix(as.numeric(as.matrix(data[1:NegContrNo,])), nrow=NegContrNo)} ## need to change here
  else {neg <- matrix(as.numeric(as.matrix(data[(PosContrNo+1):(PosContrNo+NegContrNo),])), nrow=NegContrNo)}

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
  datin <- data.frame("ecs"=as.character(c(1:length(delta_cpm))),  "uai"=delta_cpm)
  ggplot(datin, aes(x=ecs, y=uai)) +
    ylim(y_down, y_up) +
    geom_vline(xintercept=c(1:length(delta_cpm)), color = "grey85", size=0.5) +
    geom_hline(yintercept=2*stdev_cpm, linetype="dashed", color = "green", size=0.5) +
    geom_hline(yintercept=-2*stdev_cpm, linetype="dashed", color = "green", size=0.5) +
    geom_hline(yintercept=3*stdev_cpm, linetype="dashed", color = "red", size=0.5) +
    geom_hline(yintercept=-3*stdev_cpm, linetype="dashed", color = "red", size=0.5) +
    theme_bw() +
    ylab("Difference from expected value") +
    theme(axis.title.x=element_blank()) +
    scale_x_discrete(breaks=datin$ecs, labels= sample_names) +
    geom_point(color="blue4") +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  dev.off()
  #percentage of reads allocated to negative control
  perc <- colSums(neg)/lib_size
  res <- rep("OK", ncol(neg))
  res[(abs(delta_cpm) > 2*stdev_cpm) & (abs(delta_cpm) < 3*stdev_cpm)] <- "ALERT"
  res[abs(delta_cpm) >= 3*stdev_cpm] <- "FAIL"
  res_perc <- rep("OK", ncol(neg))
  res_perc[perc > 0.1] <- "ALERT"
  res_final <- rep("OK", ncol(neg))
  res_final[res == "FAIL" | res_perc =="FAIL"] <- "FAIL"
  res_final[res == "ALERT" & res_perc !="FAIL"] <- "ALERT"
  res_final[res != "FAIL" & res_perc =="ALERT"] <- "ALERT"
  table <- as.data.frame(cbind("sample_names"=sample_names, "PercentageOverLibrary"=round(perc*100, digits=3), "QC_result_deviance"=res,  "QC_result_Percentage"=res_perc, "QC_result"=res_final))
  write.table(table, "QC_neg_table.txt", quote=F, row.names=F, col.names=T)


  ## positive controls ##
  if (OrderContr == "positive") {pos <- matrix(as.numeric(as.matrix(data[1:PosContrNo,])), nrow=PosContrNo)}
  else {pos <- matrix(as.numeric(as.matrix(data[(NegContrNo+1):(NegContrNo+PosContrNo),])), nrow=PosContrNo)}

  perc2 <- colSums(pos)/lib_size
  res2 <- rep("OK", ncol(pos))
  res2[perc2 > 0.1] <- "ALERT"
  res2[perc2 > 0.4] <- "FAIL"
  table2 <- as.data.frame(cbind("sample_names"=sample_names, "PercentageOverLibrary"=round(perc2*100, digits=3), "QC_result"=res2))
  write.table(table2, "QC_pos_table.txt", quote=F, row.names=F)
  png("QC_pos_plot.png", width = 700, height = 700)

  datin2 <- data.frame("ecs"=as.character(c(1:length(delta_cpm))),  "uai"=perc2*100)
  ggplot(datin2, aes(x=ecs, y=uai)) +
    geom_vline(xintercept=c(1:length(delta_cpm)), color = "grey85", size=0.5) +
    geom_hline(yintercept=40, linetype="dashed", color = "red", size=0.5) +
    theme_bw() +
    ylab("% allocated to positive controls") +
    theme(axis.title.x=element_blank()) +
    scale_x_discrete(breaks=datin$ecs, labels= sample_names) +
    scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,100)) +
    geom_point(color="blue4") +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  dev.off()

  qc.alert <- unique(c(table$sample_names[table$QC_result == "ALERT"], table2$sample_names[table2$QC_result == "ALERT"]))
  qc.fail <- unique(c(table$sample_names[table$QC_result == "FAIL"], table2$sample_names[table2$QC_result == "FAIL"]))
  return(qc.alert)
  return(qc.fail)
}
