
normaliseHTG <- function(data.in){
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
