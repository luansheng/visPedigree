sortped <- function(ped,addgen=TRUE) {
  ped_new <- copy(ped)
  ped_new_colnames <- colnames(ped_new)
  if ("Gen" %in% ped_new_colnames) {
    ped_new[,Gen:=NULL]
  }
  sire_dam_v <- unique(c(ped_new$Sire,ped_new$Dam))
  sire_dam_v <- sire_dam_v[!is.na(sire_dam_v)]
  ped_parent_dt <- ped_new
  ind_trace_gen_dt <- setDT(list(Ind=ped_parent_dt$Ind,TraceGen=rep(0,nrow(ped_parent_dt))))
  while (sum(ped_parent_dt$Ind %chin% sire_dam_v) > 0) {
    #=== Detect pedigree loop =========================================================
    # The progeny's pedigree is subsetted
    ped_progeny_dt <- ped_parent_dt[!(Ind %chin% sire_dam_v)]
    # No offspring are subsetted because pedigree loops may cause IDs of Ind column are equal
    # to IDs of the Sire and Dam columns
    if (nrow(ped_progeny_dt) == 0) {
      stop("Pedigree error! Pedigree loops are detected!")
    }
    # Add tracing generation number
    ind_trace_gen_dt[Ind %chin% sire_dam_v,TraceGen:=TraceGen+1]
    # The parents' pedigree are subsetted
    ped_parent_dt <- ped_parent_dt[Ind %chin% sire_dam_v]
    sire_dam_v <- unique(c(ped_parent_dt$Sire,ped_parent_dt$Dam))
    sire_dam_v <- sire_dam_v[!is.na(sire_dam_v)]
    if (length(sire_dam_v[!is.na(sire_dam_v)])==0) {
      break
    }
  }

  # Assigning the progenis with parents but without progeny to the minimum tracing generation of parents - 1
  ped_trace_gen_dt <- merge(ped_new,ind_trace_gen_dt,by=c("Ind"),all.x=TRUE)
  setnames(ped_trace_gen_dt,old=c("TraceGen"),new=c("TraceGenInd"))
  TraceGenSire = TraceGenDam = NULL
  ped_trace_gen_dt[,TraceGenSire:=TraceGenInd[match(Sire,Ind)]]
  ped_trace_gen_dt[,TraceGenDam:=TraceGenInd[match(Dam,Ind)]]
  ped_trace_gen_dt[(TraceGenInd==0) & ((!is.na(Sire)) | (!is.na(Dam))),
                 TraceGenInd:=(apply(as.matrix(.SD),1,function(x) min(x,na.rm=TRUE))-1),
                 .SDcols=c("TraceGenSire","TraceGenDam")]

  # Setting the individuals witout parents and progenies as founders
  max_trace_gen_num_s <- max(ind_trace_gen_dt$TraceGen,na.rm = TRUE)
  ped_trace_gen_dt[(TraceGenInd==0) & (is.na(Sire) & is.na(Dam)),TraceGenInd:=max_trace_gen_num_s]

  # full-sib individuals have the same tracing generation number
  ped_trace_gen_dt[!is.na(Sire) | !is.na(Dam),FamilyLabel:=paste(Sire,Dam,sep="")]
  ped_trace_gen_dt[(!is.na(Sire)) | (!is.na(Dam)), MaxTraceGen:=max(TraceGenInd,na.rm=TRUE),
                                       by=c("FamilyLabel")]
  ped_trace_gen_dt[!is.na(FamilyLabel),TraceGenInd:=MaxTraceGen]

  # if an individual has not parents, it's generation number will be same with that of it's mater
  ind_no_parents_v <- ped_trace_gen_dt[is.na(Sire) & is.na(Dam), Ind]
  if (length(ind_no_parents_v) > 0) {
    sire_gen_dt <-
      unique(ped_trace_gen_dt[((Sire %chin% ind_no_parents_v) &
                                 (TraceGenSire < TraceGenDam)), .(Sire, TraceGenDam)])
    if (nrow(sire_gen_dt) > 0) {
      ped_trace_gen_dt[((Sire %chin% ind_no_parents_v) &
                          (TraceGenSire < TraceGenDam)), ":="(TraceGenSire = TraceGenDam)]
      ped_trace_gen_dt[match(sire_gen_dt$Sire, Ind), TraceGenInd := sire_gen_dt$TraceGenDam]
    }
    dam_gen_dt <-
      unique(ped_trace_gen_dt[((Dam %chin% ind_no_parents_v) &
                                 (TraceGenDam < TraceGenSire)), .(Dam, TraceGenSire)])
    if (nrow(dam_gen_dt) > 0) {
      ped_trace_gen_dt[((Dam %chin% ind_no_parents_v) &
                          (TraceGenDam < TraceGenSire)), ":="(TraceGenDam = TraceGenSire)]
      ped_trace_gen_dt[match(dam_gen_dt$Dam, Ind), TraceGenInd := dam_gen_dt$TraceGenSire]
    }
  }

  # The tracing generation number of some indivduals may be not right.
  # The following code try to renew individual tracing generation number
  # by the difference with that of parent
  # The tracing generation number of an individual will be renewed as
  # min(parents' tracing generation number)-1, if interval on tracing generation number
  # between it and it's parents is greater than 1.
  ped_trace_gen_dt[!is.na(TraceGenSire) | !is.na(TraceGenDam),
              TraceGenInterval := apply(as.matrix(.SD), 1, function(x) min(x, na.rm = TRUE)) - TraceGenInd,
              .SDcols = c("TraceGenSire", "TraceGenDam")]
  while (max(ped_trace_gen_dt$TraceGenInterval, na.rm = TRUE) > 1) {
    ped_trace_gen_dt[TraceGenInterval > 1,
                TraceGenInd := apply(as.matrix(.SD), 1, function(x) min(x, na.rm = TRUE)) - 1,
                .SDcols = c("TraceGenSire", "TraceGenDam")]
    ped_trace_gen_dt[, ":="(TraceGenSire = NULL, TraceGenDam = NULL)]
    ped_trace_gen_dt[,TraceGenSire:=TraceGenInd[match(Sire,Ind)]]
    ped_trace_gen_dt[,TraceGenDam:=TraceGenInd[match(Dam,Ind)]]
    ped_trace_gen_dt[!is.na(TraceGenSire) | !is.na(TraceGenDam),
                TraceGenInterval :=  apply(as.matrix(.SD), 1, function(x) min(x,na.rm = TRUE)) - TraceGenInd,
                .SDcols = c("TraceGenSire", "TraceGenDam")]
  }

  ped_trace_gen_dt[, ":="(TraceGenSire = NULL,
                     TraceGenDam = NULL,
                     TraceGenInterval = NULL,
                     FamilyLabel = NULL,
                     MaxTraceGen = NULL)]
  ped_trace_gen_dt[, TraceGenInd := TraceGenInd + 1]
  max_trace_gen <- max(ped_trace_gen_dt$TraceGenInd)
  # convert the tracing generation to real generation
  ped_trace_gen_dt[, Gen := (-1) * TraceGenInd + max_trace_gen + 1]
  ped_trace_gen_dt[,TraceGenInd:=NULL]
  ped_column_name <- colnames(ped_trace_gen_dt)
  ped_column_name_new <-
    c(c("Ind", "Sire", "Dam"), ped_column_name[-which(ped_column_name %chin% c("Ind", "Sire", "Dam"))])
  ped_new <-
    ped_trace_gen_dt[order(Gen,Ind), ..ped_column_name_new]
  if (!addgen) {
    ped_new[,Gen:=NULL]
  }
  return(ped_new)
}


