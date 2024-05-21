#' Read data from files of (m)any format(s)
#'
#' Wrapper for the most common `read_*()` functions from the tidyverse and
#' readxl packages that attempts to guess the file type based on the extension.
#' Currently reads `.csv` as comma separated, `.txt` as tab delimited, and
#' `.xlsx` as Excel spreadsheet. Any further arguments are passed directly to
#' the corresponding function.
#'
#' @param filename String providing the path to the input file.
#'
#' @return Data frame created by importing the file.
#' @export
#'
read_any <- function(filename, ...){
  ## Check that that the filename exists
  if(!file.exists(filename))
    stop("There is no file called ",filename,".\n")
  
  ## Extract file extension
  ext <- str_extract(filename,"\\.[0-9a-z]+$")
  
  ## Choose file reader
  myreader <- switch(ext,
                     .csv = read_csv,
                     .txt = read_tsv,
                     .xlsx = read_xlsx)
  
  ## Check if file can be read
  if(is.null(myreader))
    stop("Sorry, I don't know how to read files of type ",ext,".\n")
    
  ## Read file
  myreader(filename,...)
}