sortped <- function(ped) {
  ped_new <- copy(ped)
  sires <- unique(ped_new$Sire)
  dams  <- unique(ped_new$Dam)
  ped_new_colnames <- colnames(ped_new)
  if ("Gen" %in% ped_new_colnames) {
    ped_new[,Gen:=NULL]
  }
  sire_dam_vect <- unique(c(sires, dams))
  ped_offspring_list <- vector("list", length(sire_dam_vect))
  ped_parents <- ped_new
  i <- 1
  while (sum(ped_parents$Ind %chin% sire_dam_vect) > 0) {
    # The Offspring' pedigree are subsetted
    ped_offspring_1 <- ped_parents[!(Ind %chin% sire_dam_vect)]
    # No offspring are subsetted because pedigree loops may cause IDs of Ind column are equal
    # to IDs of the Sire and Dam columns
    if (nrow(ped_offspring_1) == 0) {
      stop("Pedigree error! Pedigree loops were detected!")
    }
    # The parents' pedigree are subsetted
    ped_parents <- ped_parents[Ind %chin% sire_dam_vect]
    # Delete individuals without parents and offspring simultaneously from the offspring pedigree
    ped_offspring_2 <-
      unique(rbind(ped_offspring_1[Sire %chin% ped_parents$Ind], ped_offspring_1[Dam %chin% ped_parents$Ind]))
    ped_offspring_list[[i]] <-
      cbind(ped_offspring_2, Gen = rep(i, nrow(ped_offspring_2)))
    sire_dam_vect <- unique(c(ped_parents$Sire, ped_parents$Dam))
    sire_dam_vect <- sire_dam_vect[!is.na(sire_dam_vect)]
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
    merge(ped_new,
          ped_new[, .(Ind, Gen)],
          by.x = "Sire",
          by.y = "Ind",
          all.x = TRUE)
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
  # The generation number of some indivduals may be not right.
  # The individuals which have no progeny were classified in the maximum generation.
  # The following code try to renew these individuals' generation number by their parent's generation number.
  # The right interval on generation number between focus individual with their parents should be equal to 1.
  # If the interval is > 1, the generation number of the foucus individual will be renewed as
  # max(parents' generation number)+1.
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
  return(ped_new)
}
