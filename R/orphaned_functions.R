#' Title
#' 
#' This function takes the number of carbons and double
# bonds and returns its respective lipid. Only works for basic
# ones
#'
#' @param carbons 
#' @param bonds 
#' @param add_ion 
#'
#' @return Character vector with entries Name and Formula.
#'
get_PC_formula <- function(carbons, bonds, add_ion = TRUE){
  C <- carbons + 8 + add_ion*2
  H <- (carbons + 8 - bonds)*2 + add_ion*3
  O <- 8 + add_ion*2
  N <- 1
  P <- 1
  name <- paste("PC ", carbons, ":", bonds, sep = "")
  formula <- paste("C", C,
                   "H", H, "N",
                   "O", O, "P", sep = "")
  lipidInfo <- c(Name=name, Formula=formula)
  return(lipidInfo)
}
