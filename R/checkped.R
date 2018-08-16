################################################
#Adapted from part of the 'prepPed' function
# written by Matthew Wolak
#in the 'nadiv' package
################################################

#' @import data.table
#' @importFrom matrixStats rowMaxs

checkped <- function(ped, sex = NULL) {
  ped_new <- copy(ped)
  ped_is_DT <- "data.table" %in% class(ped_new)
  if (!ped_is_DT) {
    setDT(ped_new)
  }
  setnames(ped_new,
           old = colnames(ped_new)[1:3],
           new = c("Ind", "Sire", "Dam"))
  setkey(ped_new, Ind, Sire, Dam)
  #===detecting and setting missing values===========================================
  ped_new[, ":="(Ind = as.character(Ind),
                 Sire = as.character(Sire),
                 Dam = as.character(Dam))]
  #Individuals will be deleted if there are "", " ", "0", "*", "NA", and NA in the Ind column.
  if (any(ped_new$Ind %in% c("", " ", "0", "*", "NA", NA))) {
    warning("There are missing values in the inividual column. These records are discarded.")
    warning("Please ensure first three columns of the pedigree are individual, sire, and dam.")
    ped_new <-
      ped_new[-which(ped_new$Ind %in% c("", " ", "0", "*", "NA", NA))]
  }
  #Missing parents are shown by "", " ", "0", "*", and "NA" in the pedigree file, they are setted as NA
  if (length(ped_new[Sire %in% c("", " ", "0", "*", "NA"), Sire]) > 0) {
    ped_new[Sire %in% c("", " ", "0", "*", "NA"), Sire := NA]
    warning("blank, Zero, asterisk, and NA means a missing parent in the sire column.")
  }

  if (length(ped_new[Dam %in% c("", " ", "0", "*", "NA"), Dam]) > 0) {
    ped_new[Dam %in% c("", " ", "0", "*", "NA"), Dam := NA]
    warning("blank, Zero, asterisk, and NA means a missing parent in the dam column.")
  }
  #The programme will stop if there are no parents in the sire and dam columns.
  if (all(is.na(ped_new$Sire)) & all(is.na(ped_new$Dam))) {
    stop("All dams and sires are missing")
  }


  #====Setting duplicated records===============================================
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
        "duplicated individual, sire, and dam are deleted in the pedigree. Only the first",
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
        "duplicated individual are found in the pedigree. Only the first",
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

  #===Checking parents' sex
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
    warning(paste(bisexual_parents, collapse = ", "))
  }

  #===renewing pedigree by adding missing parents or founders=======================
  sires_dams <- unique(c(sires, dams))
  if (sum(!(sires_dams %chin% ped_new$Ind)) > 0) {
    sires_dams_missing <- sires_dams[!(sires_dams %chin% ped_new$Ind)]
    sires_dams_missing_DT <- setDT(list(
      Ind = sires_dams_missing,
      Sire = rep(NA, length(sires_dams_missing)),
      Dam = rep(NA, length(sires_dams_missing))
    ))
    ped_new <- rbind(sires_dams_missing_DT, ped_new, fill = TRUE)
  }


  #===sorting parents in front of offspring in the individual column.
  ped_new[, SeqNum := seq(nrow(ped_new))]
  ped_new_test <-
    merge(ped_new,
          ped_new[, .(Ind, SeqNum)],
          by.x = "Sire",
          by.y = "Ind",
          all.x = TRUE)
  ped_new[, SeqNum := NULL]
  setnames(ped_new_test,
           c("SeqNum.x", "SeqNum.y"),
           c("SeqNum", "SeqNumSire"))
  ped_new_test <-
    merge(
      ped_new_test,
      ped_new_test[, .(Ind, SeqNum)],
      by.x = "Dam",
      by.y = "Ind",
      all.x = TRUE
    )
  setnames(ped_new_test,
           c("SeqNum.x", "SeqNum.y"),
           c("SeqNum", "SeqNumDam"))
  if (any(
    c(
      ped_new_test$SeqNum < ped_new_test$SeqNumSire,
      ped_new_test$SeqNum < ped_new_test$SeqNumDam
    ),
    na.rm = TRUE
  )) {
    ped_parents <- ped_new
    ped_offspring_list <- vector("list", length(sires_dams))
    i <- 1
    while (sum(ped_parents$Ind %chin% sires_dams) > 0) {
      ped_tmp_1 <- ped_parents[!(Ind %chin% sires_dams)]
      ped_parents <- ped_parents[Ind %chin% sires_dams]
      ped_tmp_2 <-
        unique(rbind(ped_tmp_1[Sire %chin% ped_parents$Ind], ped_tmp_1[Dam %chin% ped_parents$Ind]))
      ped_offspring_list[[i]] <-
        cbind(ped_tmp_2, Gen = rep(i, nrow(ped_tmp_2)))
      sires_dams <- unique(c(ped_parents$Sire, ped_parents$Dam))
      sires_dams <- sires_dams[!is.na(sires_dams)]
      i <- i + 1
    }
    ped_offspring_list[[i]] <-
      cbind(ped_parents, Gen = rep(i, nrow(ped_parents)))
    ped_offspring_DT <- do.call("rbind", ped_offspring_list)

    if (sum(!(ped_new$Ind %chin% ped_offspring_DT$Ind)) > 0) {
      ped_other <- ped_new[!(Ind %chin% ped_offspring_DT$Ind)]
      ped_other <- cbind(ped_other, Gen = rep(i, nrow(ped_other)))
      ped_new <- rbind(ped_offspring_DT, ped_other)
    } else {
      ped_new <- ped_offspring_DT
    }
    ped_new[, Gen := (-1) * Gen + i + 1]
    ped_new_Gen <-
      merge(
        ped_new,
        ped_new[, .(Ind, Gen)],
        by.x = "Sire",
        by.y = "Ind",
        all.x = TRUE
      )
    ped_new_Gen <-
      merge(
        ped_new_Gen,
        ped_new[, .(Ind, Gen)],
        by.x = "Dam",
        by.y = "Ind",
        all.x = TRUE
      )
    setnames(ped_new_Gen,
             c("Gen.x", "Gen.y", "Gen"),
             c("Gen", "GenSire", "GenDam"))
    ped_new_Gen[!is.na(GenSire) | !is.na(GenDam),
                GenInterval := Gen - rowMaxs(as.matrix(.SD), na.rm = TRUE),
                .SDcols = c("GenSire", "GenDam")]
    # The Generation number of some indivduals may be not right. The individuals which have no progeny were classified in the maximum generation. The following code try to renew these individuals' Generation number by their parent's generation number. The interval on generation number between focus individual with their parents should be equal to 1. If the interval is > 1, the generation number of the foucus individual will be renewed as max(parents' generation number)+1.
    while (max(ped_new_Gen$GenInterval, na.rm = TRUE) > 1) {
      ped_new_Gen[GenInterval > 1,
                  Gen := rowMaxs(as.matrix(.SD), na.rm = TRUE) + 1,
                  .SDcols = c("GenSire", "GenDam")]
      ped_new_Gen[, ":="(GenSire = NULL, GenDam = NULL)]
      ped_new_Gen <- merge(
        ped_new_Gen,
        ped_new_Gen[, .(Ind, Gen)],
        by.x = "Sire",
        by.y = "Ind",
        all.x = TRUE
      )
      setnames(ped_new_Gen, c("Gen.x", "Gen.y"), c("Gen", "GenSire"))
      ped_new_Gen <- merge(
        ped_new_Gen,
        ped_new_Gen[, .(Ind, Gen)],
        by.x = "Dam",
        by.y = "Ind",
        all.x = TRUE
      )
      setnames(ped_new_Gen, c("Gen.x", "Gen.y"), c("Gen", "GenDam"))
      ped_new_Gen[!is.na(GenSire) | !is.na(GenDam),
                  GenInterval := Gen - rowMaxs(as.matrix(.SD), na.rm = TRUE),
                  .SDcols = c("GenSire", "GenDam")]
    }
    ped_new_Gen[, ":="(GenSire = NULL,
                       GenDam = NULL,
                       GenInterval = NULL)]
    ped_column_name <- colnames(ped_new_Gen)
    ped_column_name_new <-
      c(c("Ind", "Sire", "Dam"), ped_column_name[-which(ped_column_name %chin% c("Ind", "Sire", "Dam"))])
    ped_new <-
      ped_new_Gen[order(Gen), ped_column_name_new, with = FALSE]
    ped_new[, Gen := NULL]
  }

  #Add each individual sex when sex parameter is not NULL
  if (!is.null(sex) & any(!is.na(ped_new$Sire))) {
    ped_new[Ind %chin% Sire,Sex:="male"]
  }
  if (!is.null(sex) & any(!is.na(ped_new$Dam))) {
    ped_new[Ind %chin% Dam,Sex:="female"]
  }
  ped_new[,Sex:=tolower(Sex)]
  return(ped_new)
}
