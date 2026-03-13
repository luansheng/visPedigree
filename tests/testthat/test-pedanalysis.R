library(testthat)
library(visPedigree)
library(data.table)

# Create a minimal testing dataset
test_ped_dt <- data.table::data.table(
  Ind = c("A", "B", "C", "D", "E", "F"),
  Sire = c(NA, NA, "A", "A", "C", "C"),
  Dam = c(NA, NA, "B", "B", "D", "D"),
  Sex = c("male", "female", "male", "female", "male", "female"),
  Gen = c(1, 1, 2, 2, 3, 3),
  Group = c("G1", "G1", "G2", "G2", "G3", "G3")
)
# C and D are full siblings (F _ offspring = 0.25). 
# E and F are offspring of full siblings.

test_ped <- suppressMessages(tidyped(test_ped_dt, reference = c("E", "F")))

test_that("pedrel calculates relation correctly and handles boundaries", {
  # N < 2 case for G1 (after subsetting by Gen, if one was dropped, or just test N=1)
  dt_n1 <- data.table::data.table(Ind = c("A", "B", "C"), Sire = c(NA, NA, "A"), Dam = c(NA, NA, "B"), Sex = c("male", "female", "male"), Gen = c(1, 1, 2))
  ped_n1 <- suppressMessages(tidyped(dt_n1))
  
  res_n1 <- suppressWarnings(pedrel(ped_n1, by = "Gen"))
  expect_equal(res_n1[Gen == 2, NUsed], 1)
  expect_true(is.na(res_n1[Gen == 2, MeanRel]))
  
  # For test_ped, group by Gen
  rel <- pedrel(test_ped, by = "Gen")
  
  # Gen 1: A, B (Unrelated) -> 0
  expect_equal(rel[Gen == 1, MeanRel], 0)
  
  # Gen 2: C, D (full sibs). a_ij = 0.5. Mean of off-diagonals = 0.5
  expect_equal(rel[Gen == 2, MeanRel], 0.5)
  
  # Test sample size behavior (less than 2 individuals after filtering)
  rel_sample <- suppressWarnings(pedrel(test_ped, by = "Gen", reference = c("A", "C")))
  expect_true(all(is.na(rel_sample$MeanRel)))
  
  # Test with a different by parameter (e.g. BirthYear)
  test_ped$YearGroup <- c(1, 1, 2, 2, 3, 3) 
  rel_year <- pedrel(test_ped, by = "YearGroup")
  expect_true("YearGroup" %in% names(rel_year))
  expect_false("Group" %in% names(rel_year))
  expect_equal(rel_year[YearGroup == 2, MeanRel], 0.5) # C, D in YearGroup 2
  
  # Test reference parameter filtering with another by variable
  # If reference is only C, D, E, then YearGroup 3 should only have E left -> NUsed=1, NA MeanRel
  rel_year_ref <- suppressWarnings(pedrel(test_ped, by = "YearGroup", reference = c("C", "D", "E")))
  expect_equal(rel_year_ref[YearGroup == 2, NUsed], 2)
  expect_equal(rel_year_ref[YearGroup == 2, MeanRel], 0.5)
  expect_equal(rel_year_ref[YearGroup == 3, NUsed], 1)
  expect_true(is.na(rel_year_ref[YearGroup == 3, MeanRel]))
})

test_that("pedrel captures deep inbreeding and correct ancestor tracing", {
  # Deep inbreeding: Full-sib mating (Gen 3), then their offspring full-sib mating (Gen 4)
  dt_deep <- data.table::data.table(
    Ind  = c("A", "B", "C", "D", "E", "F", "G", "H"),
    Sire = c(NA, NA, "A", "A", "C", "C", "E", "E"),
    Dam  = c(NA, NA, "B", "B", "D", "D", "F", "F"),
    Gen  = c(1, 1, 2, 2, 3, 3, 4, 4)
  )
  tp_deep <- tidyped(dt_deep)
  rel_deep <- pedrel(tp_deep, by = "Gen")
  
  # Gen 1: Unrelated founders
  expect_equal(rel_deep[Gen == 1, MeanRel], 0.0)
  
  # Gen 2: C, D are full sibs (0.5)
  expect_equal(rel_deep[Gen == 2, MeanRel], 0.5)
  
  # Gen 3: E, F's parents are C, D (full sibs). a_EF = 0.5 + 0.5 * f_CD = 0.5 + 0.5 * 0.5 = 0.75
  # Error check: if ancestor tracing failed, this would be 0.5
  expect_equal(rel_deep[Gen == 3, MeanRel], 0.75)
  
  # Gen 4: G, H's parents are E, F. a_GH = 0.5 + 0.5 * f_EF = 0.5 + 0.5 * 1.0 (since F_E=0.25+0.25*0.5=0.375? No, a_EF=0.75)
  # Actually a_GH = 0.5 + 0.5 * a_EF = 0.5 + 0.5 * 0.75 = 0.875? 
  # Wait, manual calibration: 
  # A_mat[G, H] = 0.5 * (A[E, G] + A[F, G]) = 0.5 * (0.5*(A[E,E]+A[E,F]) + 0.5*(A[F,E]+A[F,F]))
  # A[E,E] = 1 + f_CD = 1 + 0.25 = 1.25. A[E,F] = 0.75.
  # A_mat[G,H] = 0.5 * (0.5*(1.25+0.75) + 0.5*(0.75+1.25)) = 0.5 * (1.0 + 1.0) = 1.0
  expect_equal(rel_deep[Gen == 4, MeanRel], 1.0)
})

test_that("pedgenint computes Average from all parent-offspring pairs", {
  test_ped$BirthYear <- c(2000, 2001, 2005, 2006, 2010, 2012)
  # SS: C-A (5), E-C (5) -> mean=5

  # SD: D-A (6), F-C (7) -> mean=6.5
  # DS: C-B (4), E-D (4) -> mean=4
  # DD: D-B (5), F-D (6) -> mean=5.5
  # Average from ALL 8 pairs: (5+6+5+7+4+5+4+6)/8 = 5.25
  
  suppressMessages(
    genint_res <- pedgenint(test_ped, timevar = "BirthYear", unit = "year")
  )
  
  avg_res <- genint_res[Pathway == "Average"]
  expect_equal(avg_res$Mean, 5.25)
  expect_equal(avg_res$N, 8L)
  expect_true(!is.na(avg_res$SD))
  
  # Sex-specific pathways should still work
  expect_equal(genint_res[Pathway == "SS", Mean], 5)
  expect_equal(genint_res[Pathway == "SD", Mean], 6.5)
  expect_equal(genint_res[Pathway == "DS", Mean], 4)
  expect_equal(genint_res[Pathway == "DD", Mean], 5.5)
  
  # SO/DO: sex-independent pathways
  # SO (Sire→Offspring): A→C(5), A→D(6), C→E(5), C→F(7) -> N=4, Mean=5.75
  so_res <- genint_res[Pathway == "SO"]
  expect_equal(so_res$N, 4L)
  expect_equal(so_res$Mean, 5.75)
  
  # DO (Dam→Offspring): B→C(4), B→D(5), D→E(4), D→F(6) -> N=4, Mean=4.75
  do_res <- genint_res[Pathway == "DO"]
  expect_equal(do_res$N, 4L)
  expect_equal(do_res$Mean, 4.75)
  
  # All 7 pathways present
  expect_equal(sort(genint_res$Pathway),
               c("Average", "DD", "DO", "DS", "SD", "SO", "SS"))
})

test_that("pedcontrib f_e and f_a compute on full sets despite top cutoff", {
  cont_all <- suppressMessages(pedcontrib(test_ped, mode = "both", top = 100))
  cont_top1 <- suppressMessages(pedcontrib(test_ped, mode = "both", top = 1))
  
  # f_e and f_a should be identical because they're based on full arrays before truncation
  expect_equal(cont_all$summary$f_e, cont_top1$summary$f_e)
  expect_equal(cont_all$summary$f_a, cont_top1$summary$f_a)
  
  # Ensure reported reflects top
  expect_equal(cont_top1$summary$n_founder_show, 1)
  expect_equal(cont_all$summary$n_founder_show, length(cont_all$founders$Ind))
  expect_true(cont_top1$summary$n_founder > 1) 
})

test_that("pedpartial and pedancestry run without addnum=TRUE initially", {
  # Create tidyped without IndNum (addnum = FALSE)
  ped_nonum <- suppressMessages(tidyped(test_ped_dt, addnum = FALSE))
  
  expect_false("IndNum" %in% names(ped_nonum))
  
  # Should not crash and automatically handle missing nums
  part_res <- pedpartial(ped_nonum, ancestors = c("A", "B"), top = 2)
  expect_true(all(c("A", "B") %in% names(part_res)))
  expect_equal(nrow(part_res), nrow(ped_nonum))
  
  # pedancestry with labels
  ped_nonum$Label <- ifelse(ped_nonum$Ind %in% c("A", "B"), "Base", NA)
  anc_res <- pedancestry(ped_nonum, foundervar = "Label")
  expect_true("Base" %in% names(anc_res))
  # All non-founders descend completely from A and B, so Base should be 1.0 for all
  expect_true(all(anc_res$Base == 1.0))
})

# --- Additional tests for pedne, pedecg, pedsubpop, pedfclass ---

test_that("pedecg computes equivalent complete generations", {
  ecg_res <- pedecg(test_ped)
  expect_true(all(c("ECG", "FullGen", "MaxGen") %in% names(ecg_res)))
  expect_equal(nrow(ecg_res), nrow(test_ped))
  # Founders A, B have ECG = 0
  expect_equal(ecg_res[Ind == "A", ECG], 0)
  expect_equal(ecg_res[Ind == "B", ECG], 0)
  # C has 2 known parents -> ECG = 1
  expect_equal(ecg_res[Ind == "C", ECG], 1)
  # E has parents C, D who each have ECG=1 -> ECG = 1 + (1+1)/2 = 2
  expect_equal(ecg_res[Ind == "E", ECG], 2)
})

test_that("pedne computes effective population size by cohort", {
  test_ped$BirthYear <- c(2000, 2000, 2005, 2005, 2010, 2010)
  res <- suppressMessages(pedne(test_ped, by = "BirthYear"))
  expect_true(all(c("Cohort", "N", "DeltaC", "Ne") %in% names(res)))
  # Only cohort 2010 has ECG > 1 (founders/gen1 filtered out)
  expect_equal(nrow(res), 3)
  expect_equal(res$Cohort, c(2000, 2005, 2010))
  expect_true(any(res$DeltaC >= 0, na.rm=TRUE))
  expect_true(any(res$Ne > 0 & is.finite(res$Ne)))
})

test_that("pedsubpop splits pedigree by grouping variable", {
  res_gen <- pedsubpop(test_ped, by = "Gen")
  expect_true(is.data.table(res_gen))
  expect_equal(nrow(res_gen), length(unique(test_ped$Gen)))
  expect_true(all(c("Group", "N") %in% names(res_gen)))
  # Gen 1 (Group=1) has 2 individuals
  expect_equal(res_gen[Group == 1, N], 2)
  expect_equal(res_gen[Group == 2, N], 2)
})

test_that("pedsubpop summarizes disconnected groups and isolated individuals", {
  ped_multi_dt <- data.table::data.table(
    Ind = c("A", "B", "C", "D", "E", "F", "ISO"),
    Sire = c(NA, NA, "A", NA, NA, "D", NA),
    Dam = c(NA, NA, "B", NA, NA, "E", NA),
    Sex = c("male", "female", "male", "male", "female", "male", NA_character_)
  )
  ped_multi <- suppressMessages(tidyped(ped_multi_dt))

  res_null <- pedsubpop(ped_multi)

  expect_true(is.data.table(res_null))
  expect_equal(nrow(res_null), 3)
  expect_true(all(c("GP1", "GP2", "Isolated") %in% res_null$Group))

  expect_equal(res_null[Group == "GP1", N], 3)
  expect_equal(res_null[Group == "GP1", N_Founder], 2)

  expect_equal(res_null[Group == "GP2", N], 3)
  expect_equal(res_null[Group == "GP2", N_Founder], 2)

  expect_equal(res_null[Group == "Isolated", N], 1)
  expect_equal(res_null[Group == "Isolated", N_Sire], 0)
  expect_equal(res_null[Group == "Isolated", N_Dam], 0)
  expect_equal(res_null[Group == "Isolated", N_Founder], 1)
})

test_that("pedfclass classifies inbreeding levels", {
  res <- pedfclass(test_ped)
  expect_true(is.data.table(res))
  expect_true("FClass" %in% names(res))
  expect_s3_class(res$FClass, "ordered")
  # A, B, C, D have f=0 (4 individuals); E, F have f=0.25 (2 individuals)
  expect_equal(sum(res$Count), nrow(test_ped))
  expect_equal(res[FClass == "F = 0", Count], 4)
  expect_equal(res[FClass == "0.125 < F <= 0.25", Count], 2)
})

test_that("pedfclass handles threshold boundaries correctly", {
  threshold_ped <- data.table(
    Ind = c("A", "B", "C", "D", "E"),
    Sire = NA_character_,
    Dam = NA_character_,
    Sex = c("male", "female", "male", "female", "male"),
    f = c(0, 0.0625, 0.125, 0.25, 0.250001)
  )
  threshold_ped <- new_tidyped(threshold_ped)

  res <- pedfclass(threshold_ped)

  expect_equal(res[FClass == "F = 0", Count], 1L)
  expect_equal(res[FClass == "0 < F <= 0.0625", Count], 1L)
  expect_equal(res[FClass == "0.0625 < F <= 0.125", Count], 1L)
  expect_equal(res[FClass == "0.125 < F <= 0.25", Count], 1L)
  expect_equal(res[FClass == "F > 0.25", Count], 1L)
})

test_that("pedfclass supports custom breaks", {
  threshold_ped <- data.table(
    Ind = c("A", "B", "C", "D", "E", "F", "G"),
    Sire = NA_character_,
    Dam = NA_character_,
    Sex = c("male", "female", "male", "female", "male", "female", "male"),
    f = c(0, 0.03125, 0.0625, 0.125, 0.25, 0.5, 0.6)
  )
  threshold_ped <- new_tidyped(threshold_ped)

  res <- pedfclass(threshold_ped, breaks = c(0.03125, 0.0625, 0.125, 0.25, 0.5))

  # labels auto-generated: 5 bounded + 1 tail = 6 positives, plus "F = 0" = 7 levels
  expect_equal(
    as.character(res$FClass),
    c("F = 0", "0 < F <= 0.03125", "0.03125 < F <= 0.0625",
      "0.0625 < F <= 0.125", "0.125 < F <= 0.25", "0.25 < F <= 0.5", "F > 0.5")
  )
  expect_equal(res[FClass == "0 < F <= 0.03125",      Count], 1L)
  expect_equal(res[FClass == "0.03125 < F <= 0.0625", Count], 1L)
  expect_equal(res[FClass == "0.0625 < F <= 0.125",   Count], 1L)
  expect_equal(res[FClass == "0.125 < F <= 0.25",     Count], 1L)
  expect_equal(res[FClass == "0.25 < F <= 0.5",       Count], 1L)
  expect_equal(res[FClass == "F > 0.5",               Count], 1L)
})

test_that("pedfclass respects custom labels aligned to breaks", {
  threshold_ped <- data.table(
    Ind = c("A", "B", "C", "D", "E"),
    Sire = NA_character_,
    Dam = NA_character_,
    Sex = c("male", "female", "male", "female", "male"),
    f = c(0, 0.0625, 0.125, 0.25, 0.4)
  )
  threshold_ped <- new_tidyped(threshold_ped)

  res <- pedfclass(
    threshold_ped,
    breaks = c(0.0625, 0.125, 0.25),
    labels = c("Low", "Moderate", "High")
  )
  # tail auto-generated as "F > 0.25"
  expect_equal(as.character(res$FClass), c("F = 0", "Low", "Moderate", "High", "F > 0.25"))
  expect_equal(res[FClass == "Low",       Count], 1L)
  expect_equal(res[FClass == "High",      Count], 1L)
  expect_equal(res[FClass == "F > 0.25",  Count], 1L)

  # Wrong length should error
  expect_error(
    pedfclass(threshold_ped, breaks = c(0.0625, 0.125, 0.25), labels = c("A", "B")),
    "length equal to length\\(breaks\\)"
  )
})

test_that("pedstats returns correct structure without timevar", {
  tp <- suppressMessages(tidyped(simple_ped))
  ps <- pedstats(tp)

  # class
  expect_s3_class(ps, "pedstats")

  # $summary columns and types
  expect_true(is.data.table(ps$summary))
  expect_named(ps$summary, c("N", "NSire", "NDam", "NFounder", "MaxGen"))
  expect_equal(ps$summary$N, nrow(tp))
  expect_true(ps$summary$NFounder > 0)

  # $ecg columns
  expect_true(is.data.table(ps$ecg))
  expect_named(ps$ecg, c("Ind", "ECG", "FullGen", "MaxGen"))
  expect_equal(nrow(ps$ecg), nrow(tp))

  # no timevar -> gen_intervals is NULL
  expect_null(ps$gen_intervals)
})

test_that("pedstats returns gen_intervals with timevar", {
  tp2 <- suppressMessages(tidyped(big_family_size_ped))
  ps2 <- pedstats(tp2, timevar = "Year")

  # $gen_intervals columns
  expect_true(is.data.table(ps2$gen_intervals))
  expect_true(all(c("Pathway", "N", "Mean", "SD") %in% names(ps2$gen_intervals)))
  expect_true("Average" %in% ps2$gen_intervals$Pathway)

  # ecg = FALSE suppresses ecg
  ps_no_ecg <- pedstats(tp2, timevar = "Year", ecg = FALSE)
  expect_null(ps_no_ecg$ecg)

  # genint = FALSE suppresses gen_intervals
  ps_no_gi <- pedstats(tp2, timevar = "Year", genint = FALSE)
  expect_null(ps_no_gi$gen_intervals)
})
