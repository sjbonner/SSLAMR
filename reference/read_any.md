# Read data from files of (m)any format(s)

Wrapper for the most common `read_*()` functions from the tidyverse and
readxl packages that attempts to guess the file type based on the
extension. Currently reads `.csv` as comma separated, `.txt` as tab
delimited, and `.xlsx` as Excel spreadsheet. Any further arguments are
passed directly to the corresponding function.

## Usage

``` r
read_any(filename, show_col_types = FALSE, ...)
```

## Arguments

- filename:

  String providing the path to the input file. (string)

- show_col_types:

  If FALSE, do not show the guessed column types. (boolean)

- ...:

## Value

Data frame created by importing the file.
