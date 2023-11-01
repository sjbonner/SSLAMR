read_any <- function(filename){
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
  myreader(filename)
}