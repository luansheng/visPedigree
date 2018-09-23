#' @import data.table


# checkped will add a new column: Sex
checkped <- function(ped,addgen=TRUE) {
  ped_new <- copy(ped)
  ped_is_DT <- "data.table" %in% class(ped_new)
  if (!ped_is_DT) {
    setDT(ped_new)
  }
  setnames(ped_new,
           old = colnames(ped_new)[1:3],
           new = c("Ind", "Sire", "Dam"))
  setkey(ped_new, Ind, Sire, Dam)
  #===detecting and setting missing values=============================================
  Ind <- Sire <- Dam <- NULL
  ped_new[, ":="(Ind = as.character(Ind),
                 Sire = as.character(Sire),
                 Dam = as.character(Dam))]
  #Individuals will be deleted if there are "", " ", "0", "*", "NA", and NA in the Ind column.
  if (any(ped_new$Ind %in% c("", " ", "0", "*", "NA", NA))) {
    warning("There are missing values in the inividual column of the pedigree. These records are discarded.")
    warning("Individual, sire, and dam must be the first three columns of the pedigree.")
    ped_new <-
      ped_new[-which(ped_new$Ind %in% c("", " ", "0", "*", "NA", NA))]
  }
  # Missing parents are shown by "", " ", "0", "*", and "NA" in the pedigree file,
  # they are setted as NA
  if (length(ped_new[Sire %in% c("", " ", "0", "*", "NA"), Sire]) > 0) {
    ped_new[Sire %in% c("", " ", "0", "*", "NA"), Sire := NA]
    warning("Blank, Zero, asterisk, or NA are read as a missing parent in the sire column of the pedigree.")
  }

  if (length(ped_new[Dam %in% c("", " ", "0", "*", "NA"), Dam]) > 0) {
    ped_new[Dam %in% c("", " ", "0", "*", "NA"), Dam := NA]
    warning("Blank, Zero, asterisk, or NA are read as a missing parent in the dam column of the pedigree.")
  }
  #The programme will stop if there are no parents in the sire and dam columns.
  if (all(is.na(ped_new$Sire)) & all(is.na(ped_new$Dam))) {
    stop("All dams and sires are missing! No pedigee! Please check it!")
  }


  #====detect or delete duplicated records===============================================
  # If the duplicated records have the same individual, sire, and dam ID,
  # one record will be kept, other records will be deleted.
  if (anyDuplicated(ped_new, by = c("Ind", "Sire", "Dam")) > 0) {
    ped_new_dup <-
      ped_new[duplicated(ped_new, by = c("Ind", "Sire", "Dam"))]
    ped_new_dup_num <- nrow(ped_new_dup)
    n = ped_new_dup_num
    if (n > 5) {
      n = 5
    }
    warning(
      paste(
        "The",
        ped_new_dup_num,
        "duplicated individual, sire, and dam IDs are deleted in the pedigree. Only the first",
        n,
        "records are shown"
      )
    )
    for (i in 1:n) {
      warning(paste(ped_new_dup[i], collapse = ", "))
    }
    if (ped_new_dup_num > 5) {
      warning("...")
    }

    ped_new <- unique(ped_new, by = c("Ind", "Sire", "Dam"))
  }
  # If the duplicated records only have the same individual ID,
  #  and their sire and dam ID are different,
  # this program will stop because it is a fatal error.
  if (anyDuplicated(ped_new, by = c("Ind")) > 0) {
    ped_new_dup <-
      ped_new[duplicated(ped_new, by = c("Ind"))]
    ped_new_dup_num <- nrow(ped_new_dup)
    n = ped_new_dup_num
    if (n > 5) {
      n = 5
    }
    warning(
      paste(
        "The",
        ped_new_dup_num,
        "duplicated individual IDs are found in the pedigree. Only the first",
        n,
        "records are shown."
      )
    )

    for (i in 1:n) {
      warning(paste(ped_new_dup[i], collapse = ", "))
    }
    if (ped_new_dup_num > 5) {
      warning("...")
    }
    stop("Please check the pedigree!")
  }

  #===find bisexual parents============================================================
  sires <- unique(ped_new$Sire)
  if (any(is.na(sires))) {
    sires <- sires[-which(is.na(sires))]
  }

  dams  <- unique(ped_new$Dam)
  if (any(is.na(dams))) {
    dams  <- dams[-which(is.na(dams))]
  }
  bisexual_parents <- sires[sires %chin% dams]
  if (length(bisexual_parents) > 0) {
    warning("The following individuals are simultaneously bisexual.")
    warning(paste(sort(unique(bisexual_parents)), collapse = ", "))
  }

  #===renewing pedigree by adding missing parents or founders==========================
  sire_dam_vect <- unique(c(sires, dams))
  if (sum(!(sire_dam_vect %chin% ped_new$Ind)) > 0) {
    sire_dam_vect_missing <- sire_dam_vect[!(sire_dam_vect %chin% ped_new$Ind)]
    sire_dam_vect_missing_DT <- setDT(list(
      Ind = sire_dam_vect_missing,
      Sire = rep(NA, length(sire_dam_vect_missing)),
      Dam = rep(NA, length(sire_dam_vect_missing))
    ))
    ped_new <- rbind(sire_dam_vect_missing_DT, ped_new, fill = TRUE)
  }

  #===sorting parents in front of offspring in the individual column.
  SeqNumInd <- SeqNumSire <- SeqNumDam <- NULL
  ped_new[, SeqNumInd := .I]
  ped_new[,SeqNumSire:=SeqNumInd[match(Sire,Ind)]]
  ped_new[,SeqNumDam:=SeqNumInd[match(Dam,Ind)]]
  # Individuals are resorted  if their order are not right.
  if (any(
    c(
      ped_new$SeqNumInd < ped_new$SeqNumSire,
      ped_new$SeqNumInd < ped_new$SeqNumDam
    ),
    na.rm = TRUE
  ) | addgen)  {
    ped_new <- sortped(ped_new,addgen)
  }

  # delete internal fields
  ped_new[,":="(SeqNumInd=NULL,SeqNumSire=NULL,SeqNumDam=NULL)]

  #===Add individual sex==========================================================
  col_names <- colnames(ped_new)
  Sex <- NULL
  if (!("Sex" %in% col_names)) {
    if (any(!is.na(ped_new$Sire))) {
      ped_new[Ind %chin% Sire, Sex:="male"]
    }
    if (any(!is.na(ped_new$Dam))) {
      ped_new[Ind %chin% Dam, Sex:="female"]
    }
  }
  if ("Sex" %in% col_names) {
    if (length(ped_new[Sex %in% c("", " ", "NA"), Sex]) > 0) {
      ped_new[Sex %in% c("", " ", "NA"), Sex := NA]
      warning("Blank and NA are recoded as a missing sex in the Sex column of the pedigree.")
    }

    if (any(!is.na(ped_new$Sire))) {
      ped_new[is.na(Sex) & (Ind %chin% Sire), Sex:="male"]
    }
    if (any(!is.na(ped_new$Dam))) {
      ped_new[is.na(Sex) & (Ind %chin% Dam), Sex:="female"]
    }
  }
  if (any(c(!is.na(ped_new$Sire),!is.na(ped_new$Dam)))) {
    ped_new[!is.na(Sex),Sex:=tolower(Sex)]
    sex_name <- unique(ped_new[!is.na(Sex),Sex])
    if (!all(sex_name %in% c("male","female"))) {
      message("There are other sexes except male and female in the Sex column!")
      message(paste("Sex:",paste(sex_name,collapse = " "),sep=" "))
    }
  }
  return(ped_new)
}
