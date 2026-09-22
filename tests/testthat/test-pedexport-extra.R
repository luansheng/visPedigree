# Supplementary systematic tests for pedexport(): format invariants across
# datasets, file round-trips, xref consistency, and custom missing symbols.
# Complements test-pedexport.R (which covers API behaviour on small examples).

library(data.table)

# ---- helpers -----------------------------------------------------------------

# Invariants that every numeric (integer-coded) export must satisfy.
expect_valid_numeric_export <- function(out, tp) {
  expect_identical(names(out), c("IndNum", "SireNum", "DamNum"))
  expect_true(is.integer(out$IndNum) || is.numeric(out$IndNum))
  expect_identical(out$IndNum, seq_len(nrow(tp)))
  # every known parent has a strictly smaller code than the offspring
  expect_true(all(out$SireNum == 0L | out$SireNum < out$IndNum))
  expect_true(all(out$DamNum == 0L | out$DamNum < out$IndNum))
  # xref maps back to the original IDs bijectively
  xr <- attr(out, "xref")
  expect_equal(nrow(xr), nrow(tp))
  expect_setequal(xr$Ind, tp$Ind)
  expect_identical(xr$IndNum, out$IndNum)
  # numeric codes identify the same parents as the character IDs
  sire_char <- xr$Ind[match(out$SireNum, xr$IndNum)]
  dam_char <- xr$Ind[match(out$DamNum, xr$IndNum)]
  tp_ord <- tp[order(IndNum)]
  expect_identical(sire_char, tp_ord$Sire)
  expect_identical(dam_char, tp_ord$Dam)
}

# Invariants for character exports: parents precede offspring, founders first,
# every parent either appears as a record or is the missing symbol.
expect_valid_char_export <- function(out, tp, miss = "0") {
  expect_equal(ncol(out), 3)
  ids <- out[[1]]
  expect_setequal(ids, tp$Ind)
  sire_pos <- match(out[[2]], ids)
  dam_pos <- match(out[[3]], ids)
  is_miss_s <- if (is.na(miss)) is.na(out[[2]]) else out[[2]] == miss
  is_miss_d <- if (is.na(miss)) is.na(out[[3]]) else out[[3]] == miss
  expect_true(all(is_miss_s | sire_pos < seq_len(nrow(out))))
  expect_true(all(is_miss_d | dam_pos < seq_len(nrow(out))))
  # structure matches the source pedigree row for row
  src <- tp[match(ids, tp$Ind)]
  got_s <- ifelse(is_miss_s, NA_character_, out[[2]])
  got_d <- ifelse(is_miss_d, NA_character_, out[[3]])
  expect_identical(got_s, src$Sire)
  expect_identical(got_d, src$Dam)
}

datasets <- list(
  simple_ped = simple_ped,
  small_ped = small_ped,
  deep_ped = deep_ped,
  half_founder_ped = half_founder_ped,
  inbred_ped = inbred_ped
)

# ---- 1. invariants across datasets and formats -------------------------------

test_that("numeric formats satisfy integer-coding invariants on all datasets", {
  for (nm in names(datasets)) {
    tp <- tidyped(datasets[[nm]])
    for (sw in c("blupf90", "wombat", "mtdfreml", "dmu", "numeric")) {
      expect_valid_numeric_export(pedexport(tp, software = sw), tp)
    }
  }
})

test_that("character formats satisfy ordering invariants on all datasets", {
  for (nm in names(datasets)) {
    tp <- tidyped(datasets[[nm]])
    for (sw in c("asreml", "echidna", "hiblup")) {
      expect_valid_char_export(pedexport(tp, software = sw), tp, miss = "0")
    }
    expect_valid_char_export(pedexport(tp, software = "sommer"), tp, miss = NA)
  }
})

test_that("invariants hold for traced, bottom-aligned and inbreed pedigrees", {
  tp1 <- tidyped(simple_ped, cand = c("J5X804", "J3Y620"), trace = "all")
  tp2 <- tidyped(simple_ped, genmethod = "bottom")
  tp3 <- tidyped(inbred_ped, inbreed = TRUE)
  for (tp in list(tp1, tp2, tp3)) {
    expect_valid_numeric_export(pedexport(tp, "blupf90"), tp)
    expect_valid_char_export(pedexport(tp, "asreml"), tp, miss = "0")
  }
})

# ---- 2. file round-trips ------------------------------------------------------

test_that("numeric file round-trips through fread and tidyped identically", {
  tp <- tidyped(simple_ped)
  f <- tempfile(fileext = ".txt")
  on.exit(unlink(c(f, paste0(f, ".xref"))))
  pedexport(tp, software = "blupf90", file = f)

  re <- fread(f, header = FALSE)
  expect_equal(ncol(re), 3L)
  expect_true(all(re$V1 == seq_len(nrow(tp))))
  tp2 <- tidyped(re)
  orig <- tp[order(IndNum)]
  expect_equal(nrow(tp2), nrow(tp))
  expect_identical(ifelse(is.na(tp2$Sire), 0L, as.integer(tp2$Sire)),
                   orig$SireNum)
  expect_identical(ifelse(is.na(tp2$Dam), 0L, as.integer(tp2$Dam)),
                   orig$DamNum)

  # .xref file maps numeric codes back to original IDs
  xr <- fread(paste0(f, ".xref"), header = FALSE)
  expect_identical(xr$V1, seq_len(nrow(tp)))
  expect_setequal(xr$V2, tp$Ind)
})

test_that("character file round-trips and preserves the pedigree structure", {
  tp <- tidyped(simple_ped)
  f <- tempfile(fileext = ".txt")
  on.exit(unlink(f))
  pedexport(tp, software = "hiblup", file = f)  # no header, "0" = missing

  re <- fread(f, header = FALSE)
  expect_equal(ncol(re), 3L)
  tp2 <- tidyped(re)
  key <- c("Ind", "Sire", "Dam", "Gen")
  expect_identical(tp2[order(Ind), ..key], tp[order(Ind), ..key])
})

test_that("asreml file with header needs (and works with) an explicit skip", {
  tp <- tidyped(simple_ped)
  f <- tempfile(fileext = ".txt")
  on.exit(unlink(f))
  pedexport(tp, software = "asreml", file = f)  # header written by default

  lines <- readLines(f)
  expect_match(lines[1], "animal\\s+sire\\s+dam")
  re <- fread(f, skip = 1, header = FALSE)
  expect_equal(nrow(re), nrow(tp))
  tp2 <- tidyped(re)
  expect_setequal(tp2$Ind, tp$Ind)
})

# ---- 3. custom missing symbols ------------------------------------------------

test_that("numeric formats honour a custom integer missing symbol", {
  tp <- tidyped(simple_ped)
  out <- pedexport(tp, software = "dmu", missing = 9L)
  founders <- tp[is.na(Sire) & is.na(Dam), IndNum]
  expect_true(all(out$SireNum[founders] == 9L))
  expect_true(all(out$DamNum[founders] == 9L))
  expect_true(all(out$SireNum[-founders] != 9L | out$SireNum[-founders] == 9L &
                    out$SireNum[-founders] < out$IndNum[-founders]))
})

test_that("character formats honour a custom character missing symbol", {
  tp <- tidyped(simple_ped)
  out <- pedexport(tp, software = "asreml", missing = "NA")
  founders <- tp[is.na(Sire) & is.na(Dam), Ind]
  got <- out[match(founders, out[[1]])]
  expect_true(all(got[[2]] == "NA" & got[[3]] == "NA"))
  expect_valid_char_export(out, tp, miss = "NA")
})

# ---- 4. reconstruction paths ---------------------------------------------------

test_that("addnum/addgen FALSE reconstruction stays valid (order may differ)", {
  full <- tidyped(simple_ped)
  bare <- tidyped(simple_ped, addnum = FALSE, addgen = FALSE)

  # reconstruction produces valid exports, but the row order (and hence the
  # numeric coding) follows the topological row order of the input, which is
  # NOT the (Gen, Ind) order used when the columns are present
  num_bare <- pedexport(bare, "blupf90")
  expect_identical(num_bare$IndNum, seq_len(nrow(bare)))
  expect_true(all(num_bare$SireNum == 0L | num_bare$SireNum < num_bare$IndNum))
  expect_true(all(num_bare$DamNum == 0L | num_bare$DamNum < num_bare$IndNum))
  xr <- attr(num_bare, "xref")
  src <- bare[match(xr$Ind, bare$Ind)]
  expect_identical(xr$Ind[match(num_bare$SireNum, xr$IndNum)] , src$Sire)
  expect_identical(xr$Ind[match(num_bare$DamNum, xr$IndNum)], src$Dam)

  chr_bare <- pedexport(bare, "asreml")
  expect_valid_char_export(chr_bare, bare, miss = "0")

  # caller's object is not modified by reconstruction
  expect_false(any(c("IndNum", "Gen") %in% names(bare)))
})

# ---- 5. scale -----------------------------------------------------------------

test_that("pedexport stays correct on a 50k simulated pedigree", {
  set.seed(7)
  n_founders <- 500
  n_gen <- 8
  per_gen <- 6000
  inds <- paste0("F", seq_len(n_founders))
  sex <- rep(c("male", "female"), length.out = n_founders)
  recs <- list(data.table(Ind = inds, Sire = NA_character_, Dam = NA_character_))
  males <- inds[sex == "male"]; females <- inds[sex == "female"]
  for (g in seq_len(n_gen)) {
    new <- paste0("G", g, "_", seq_len(per_gen))
    recs[[g + 1L]] <- data.table(Ind = new,
                                 Sire = sample(males, per_gen, replace = TRUE),
                                 Dam = sample(females, per_gen, replace = TRUE))
    sx <- sample(c("male", "female"), per_gen, replace = TRUE)
    males <- c(males, new[sx == "male"]); females <- c(females, new[sx == "female"])
  }
  tp <- tidyped(rbindlist(recs)[sample(.N)])

  out <- pedexport(tp, "blupf90")
  expect_valid_numeric_export(out, tp)
  outc <- pedexport(tp, "asreml")
  expect_valid_char_export(outc, tp, miss = "0")
})
