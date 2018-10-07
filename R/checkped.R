#' Check a pedigree
#'
#' \code{checkped} function checks a pedigree.
#'
#' This function takes a pedigree, detects missaing parents, checks duplicated, bisexual individuals,  adds missing founders, and sorts the pedigree. All individuals's sex will be inferred if there is not sexual information in the pedigee. If the pedigree includes the column \strong{Sex}, then individuals's sexes need to be recoded as "male", "female", or NA (unknown sex). Missing sexes will be identified from the pedigree structure and be added if possible.
#' @param ped A data.table or data frame including the pedigree, which including the first three columns: \strong{individual}, \strong{sire} and \strong{dam} IDs. More columns, such as sex, generation can be included in the pedigree file. Names of the three columns can be assigend as you would like, but their orders must be not changded in the pedigree. Individual ID should not be coded as "", " ", "0", asterisk, and "NA", otherwise these individuals will be deleted from the pedigree. Missing parents should be denoted by either "NA", "0", asterisk. Space and "" will also be recoded as missing parents, but not be recommended.
#' @param addgen A logical value indicates whether individual generation number will be generated. The default values is TRUE, then a new column named \strong{Gen} will be added in the returned data.table.
#'
#' @return A data.table including the checked pedigree is returned. Individual, sire and dam ID columns are renamed as \strong{Ind}, \strong{Sire} and \strong{Dam}. Missing parent is replacted with the default missing value \strong{NA}. The column \strong{Sex} includes individuals' sex (male or female, NA for unknown sex). The column \strong{Gen} will be included when the parameter \emph{addgen} is TRUE. Ind, Sire, Dam and Sex columns are character; The Gen column is integer.
#' @keywords internal
#' @import data.table
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
    warning("Blank, Zero, asterisk, or character NA is recognized as a missing dam and is replaced with missing value NA in the sire column of the pedigree.")
  }

  if (length(ped_new[Dam %in% c("", " ", "0", "*", "NA"), Dam]) > 0) {
    ped_new[Dam %in% c("", " ", "0", "*", "NA"), Dam := NA]
    warning("Blank, Zero, asterisk, or NA is recognized as a missing sire and is replaced with missing value NA in the dam column of the pedigree.")
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
