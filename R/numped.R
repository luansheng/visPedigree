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
  setkey(ped_new,Ind,Sire,Dam)
  ped_new_num <- ped_new
  #ped_new_num <- merge(ped_new,ped_new[,.(Ind,IndNum)],by.x="Sire",by.y="Ind",all.x=TRUE)
  ped_new[,SireNum:=IndNum[match(Sire,Ind)]]
  #setnames(ped_new_num,c("IndNum.x","IndNum.y"),c("IndNum","SireNum"))
  #ped_new_num <- merge(ped_new_num,ped_new_num[,.(Ind,IndNum)],by.x="Dam",by.y="Ind",all.x=TRUE)
  #setnames(ped_new_num,c("IndNum.x","IndNum.y"),c("IndNum","DamNum"))
  ped_new[,DamNum:=IndNum[match(Dam,Ind)]]
  ped_new[is.na(SireNum),SireNum:=0]
  ped_new[is.na(DamNum),DamNum:=0]
  ped_column_name <- colnames(ped_new)
  ped_column_name_new <- c(c("Ind","Sire","Dam"),ped_column_name[-which(ped_column_name %chin% c("Ind","Sire","Dam"))])
  ped_new <- ped_new[, ..ped_column_name_new]
  return(ped_new[order(IndNum)])
}
