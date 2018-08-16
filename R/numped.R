numped <- function(ped) {
  ped_new <- copy(ped)
  ped_is_DT <- "data.table" %in% class(ped_new)
  if (!ped_is_DT) {
    setDT(ped_new)
  }
  setnames(ped_new,
           old = colnames(ped_new)[1:3],
           new = c("Ind", "Sire", "Dam"))
  setkey(ped_new,Ind,Sire,Dam)
  ped_new[,IndNum:=seq(nrow(ped_new))]
  ped_new_num <- merge(ped_new,ped_new[,.(Ind,IndNum)],by.x="Sire",by.y="Ind",all.x=TRUE)
  setnames(ped_new_num,c("IndNum.x","IndNum.y"),c("IndNum","SireNum"))
  ped_new_num <- merge(ped_new_num,ped_new_num[,.(Ind,IndNum)],by.x="Dam",by.y="Ind",all.x=TRUE)
  setnames(ped_new_num,c("IndNum.x","IndNum.y"),c("IndNum","DamNum"))
  ped_new_num[is.na(SireNum),SireNum:=0]
  ped_new_num[is.na(DamNum),DamNum:=0]
  ped_column_name <- colnames(ped_new_num)
  ped_column_name_new <- c(c("Ind","Sire","Dam"),ped_column_name[-which(ped_column_name %chin% c("Ind","Sire","Dam"))])
  ped_new <- ped_new_num[,ped_column_name_new,with=FALSE]
  return(ped_new[order(IndNum)])
}
