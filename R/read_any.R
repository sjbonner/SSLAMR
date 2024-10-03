#' Read data from files of (m)any format(s)
#'
#' Wrapper for the most common `read_*()` functions from the tidyverse and
#' readxl packages that attempts to guess the file type based on the extension.
#' Currently reads `.csv` as comma separated, `.txt` as tab delimited, and
#' `.xlsx` as Excel spreadsheet. Any further arguments are passed directly to
#' the corresponding function.
#'
#' @param filename String providing the path to the input file. (string)
#' @param show_col_types If FALSE, do not show the guessed column types. (boolean) 
#' @param ... 
#'
#' @return Data frame created by importing the file.
#' @export
#'
read_any <- function(filename, show_col_types = FALSE, ...){
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
  if(ext != ".xlsx")
    myreader(filename, show_col_types = show_col_types,...)
  else
    myreader(filename, ...)
}