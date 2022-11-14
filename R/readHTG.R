#'Extract count data from HTG EdgeSeq software's output file
#'
#'Prepares raw output from the HTG EdgeSeq analysis software for the quality check
#'This function takes the unformatted excel file and converts it into a dataframe to be used as input for the qualityCheck() function
#'
#'@param data.in
#'The absolute path to the excel file with HTG EdgeSeq analysis platform data
#'
#'@param PreFormatted
#'Whether the original file was pre formatted.
#'If FALSE, the excel is the output from the HTG EdgeSeq analysis platform. This includes a "Data" sheet, a box termed "Sample Names", and the fist positive control should be termed ANT1.
#'If TRUE, the preformatted excel should contain sample names in the first row, while the first column should corrspond to gene names.
#'
#'@details
#'The function formats and prepares an input excel file to be used for other functions in the package, mainly qualityCheck()
#'
#'@return
#'A dataframe to be used for other functions in the package. This still includes the positive and negative controls.
#'
#'@export
#'

readHTG <- function(data.in, PreFormatted=FALSE){

  if (PreFormatted == TRUE){

    data <- as.data.frame(read_excel(data.in))
    row.names(data) <- data[,1]
    data <- data[,-1]
  }

  else {

    # checks
    assertthat::assert_that(!is.na(match("Data", excel_sheets(data.in))), msg = "Excel file should contain a 'Data' sheet")

    # read data
    data <- as.data.frame(read_excel(data.in, "Data"))

    assertthat::assert_that(!is.na(match("Sample Name", data[,1])), msg = "File should contain a 'Sample Name' box")
    assertthat::assert_that(!is.na(match("ANT1", data[,1])), msg = "First positive control should be termed 'ANT1'")

    # get index of sample names and gene starting position
    ind.id <- match("Sample Name", data[,1])
    ind.genes <- match("ANT1", data[,1])

    # gene-only data
    samples.id <- data[ind.id,-1]
    data <- data[ind.genes:nrow(data),]
    row.names(data) <- data[,1]
    data <- data[,-1]
    data <- as.data.frame(t(apply(data, 1, function(x) as.numeric(as.character(x)))))
    colnames(data) <- samples.id
  }

  return(data)
}
