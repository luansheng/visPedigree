#' Number a pedigree
#'
#' \code{numped} function replaces individual IDs with integers in the pedigree.
#'
#' This function takes a pedigree, assigns individual IDs to continuous integers 1...n.
#' @param ped A data.table or data frame including the pedigree, in which the names of the first three columns are \strong{Ind}, \strong{Sire} and \strong{Dam}. Missing parent has been replaced with the default missing value \strong{NA}.
#'
#' @return A data.table including the new integer columns (IndNum, SireNum and DamNum) in the pedigree.
#'
#' @import data.table
numped <- function(ped) {
  ped_new <- copy(ped)
  ped_is_DT <- "data.table" %in% class(ped_new)
  if (!ped_is_DT) {
    setDT(ped_new)
  }
  setnames(ped_new,
           old = colnames(ped_new)[1:3],
           new = c("Ind", "Sire", "Dam"))
  ped_new_colnames <- colnames(ped_new)
  IndNum <- SireNum <- DamNum <- NULL
  if ("IndNum" %in% ped_new_colnames) {
    ped_new[,IndNum:=NULL]
  }
  if ("SireNum" %in% ped_new_colnames) {
    ped_new[,SireNum:=NULL]
  }
  if ("DamNum" %in% ped_new_colnames) {
    ped_new[,DamNum:=NULL]
  }
  ped_new[,IndNum:=.I]
  ped_new[,SireNum:=IndNum[match(Sire,Ind)]]
  ped_new[,DamNum:=IndNum[match(Dam,Ind)]]
  ped_new[is.na(SireNum),SireNum:=0]
  ped_new[is.na(DamNum),DamNum:=0]
  ped_column_name <- colnames(ped_new)
  ped_column_name_new <- c(c("Ind","Sire","Dam"),
                           ped_column_name[-which(ped_column_name %chin% c("Ind","Sire","Dam"))])
  ped_new <- ped_new[, ..ped_column_name_new]
  return(ped_new[order(IndNum)])
}
