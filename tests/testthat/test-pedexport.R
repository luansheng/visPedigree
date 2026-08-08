# Tests for pedexport()

library(testthat)
library(data.table)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_simple_ped <- function() {
  # A -> C, B -> C; C and E -> D (simple 3-generation chain).
  # D has distinct parents (C x E): an individual must not appear as both
  # sire and dam unless tidyped(selfing = TRUE).
  data.frame(
    Ind  = c("A", "B", "C", "D", "E"),
    Sire = c(NA,  NA,  "A", "C", NA),
    Dam  = c(NA,  NA,  "B", "E", NA),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 1. Basic structure tests
# ---------------------------------------------------------------------------

test_that("pedexport returns invisible data.table for blupf90", {
  tp  <- tidyped(make_simple_ped())
  out <- pedexport(tp, software = "blupf90")

  expect_s3_class(out, "data.table")
  expect_equal(ncol(out), 3L)
  expect_named(out, c("IndNum", "SireNum", "DamNum"))
  expect_equal(nrow(out), nrow(tp))
})

test_that("pedexport returns 3-column data.table for asreml", {
  tp  <- tidyped(make_simple_ped())
  out <- pedexport(tp, software = "asreml")

  expect_s3_class(out, "data.table")
  expect_equal(ncol(out), 3L)
  expect_named(out, c("animal", "sire", "dam"))
  expect_equal(nrow(out), nrow(tp))
})

test_that("pedexport wombat returns character columns with '0' missing", {
  tp   <- tidyped(make_simple_ped())
  womb <- pedexport(tp, software = "wombat")

  # Wombat accepts alphanumeric IDs (recoding them internally), so the
  # export keeps character IDs rather than renumbering.
  expect_named(womb, c("animal", "sire", "dam"))
  expect_type(womb$animal, "character")
  expect_type(womb$sire,  "character")

  founders <- tp[is.na(Sire) & is.na(Dam), Ind]
  expect_true(all(womb[animal %in% founders, sire] == "0"))
  expect_true(all(womb[animal %in% founders, dam]  == "0"))

  # Non-missing parents retain their character IDs
  expect_equal(womb[animal == "C", sire], "A")
  expect_equal(womb[animal == "C", dam],  "B")
})

test_that("pedexport mtdfreml is identical to blupf90", {
  tp   <- tidyped(make_simple_ped())
  b90  <- pedexport(tp, software = "blupf90")
  mtd  <- pedexport(tp, software = "mtdfreml")
  expect_equal(b90, mtd)
})

test_that("pedexport numeric returns integer columns with header by default", {
  tp  <- tidyped(make_simple_ped())
  out <- pedexport(tp, software = "numeric")
  expect_named(out, c("IndNum", "SireNum", "DamNum"))
  expect_type(out$IndNum, "integer")
})

# ---------------------------------------------------------------------------
# 2. Missing parent encoding
# ---------------------------------------------------------------------------

test_that("blupf90: founders have SireNum = 0 and DamNum = 0", {
  tp  <- tidyped(make_simple_ped())
  out <- pedexport(tp, software = "blupf90")

  founders <- tp[is.na(Sire) & is.na(Dam), IndNum]
  expect_true(all(out[IndNum %in% founders, SireNum] == 0L))
  expect_true(all(out[IndNum %in% founders, DamNum]  == 0L))
})

test_that("asreml: founders have sire = '0' and dam = '0'", {
  tp  <- tidyped(make_simple_ped())
  out <- pedexport(tp, software = "asreml")

  founders <- tp[is.na(Sire) & is.na(Dam), Ind]
  expect_true(all(out[animal %in% founders, sire] == "0"))
  expect_true(all(out[animal %in% founders, dam]  == "0"))
})

test_that("asreml: non-missing parents retain their character IDs", {
  tp  <- tidyped(make_simple_ped())
  out <- pedexport(tp, software = "asreml")

  # Individual C has parents A and B
  row_C <- out[animal == "C"]
  expect_equal(row_C$sire, "A")
  expect_equal(row_C$dam,  "B")
})

# ---------------------------------------------------------------------------
# 3. Row ordering
# ---------------------------------------------------------------------------

test_that("blupf90: rows are sorted by IndNum (parents before offspring)", {
  tp  <- tidyped(make_simple_ped())
  out <- pedexport(tp, software = "blupf90")

  expect_equal(out$IndNum, sort(out$IndNum))

  # Every non-zero SireNum must be less than the offspring IndNum
  non_founder <- out[SireNum > 0L]
  expect_true(all(non_founder$SireNum < non_founder$IndNum))
  non_founder_d <- out[DamNum > 0L]
  expect_true(all(non_founder_d$DamNum < non_founder_d$IndNum))
})

test_that("asreml: founders appear in first rows (Gen ordering)", {
  tp  <- tidyped(make_simple_ped())
  out <- pedexport(tp, software = "asreml")

  founders <- tp[is.na(Sire) & is.na(Dam), Ind]
  first_rows <- head(out$animal, length(founders))
  expect_true(all(founders %in% first_rows))
})

# ---------------------------------------------------------------------------
# 4. File output
# ---------------------------------------------------------------------------

test_that("pedexport writes a readable file for blupf90", {
  tp  <- tidyped(make_simple_ped())
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)

  pedexport(tp, software = "blupf90", file = tmp)

  expect_true(file.exists(tmp))
  lines <- readLines(tmp)
  # No header by default for blupf90
  expect_equal(length(lines), nrow(tp))
  # Each line has 3 space-separated integer tokens
  tokens <- strsplit(lines[[1]], " ")[[1]]
  expect_equal(length(tokens), 3L)
  expect_true(all(!is.na(suppressWarnings(as.integer(tokens)))))
})

test_that("asreml file writes a header by default", {
  tp  <- tidyped(make_simple_ped())
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)

  pedexport(tp, software = "asreml", file = tmp)

  lines <- readLines(tmp)
  # A header is written by default for readability; when reading the file
  # back, ASReml needs the !SKIP 1 qualifier.
  expect_equal(lines[[1]], "animal sire dam")
  expect_equal(length(lines), nrow(tp) + 1L)
})

test_that("asreml file omits the header when header = FALSE", {
  tp  <- tidyped(make_simple_ped())
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)

  pedexport(tp, software = "asreml", file = tmp, header = FALSE)

  lines <- readLines(tmp)
  expect_equal(length(lines), nrow(tp))
  expect_false(lines[[1]] == "animal sire dam")
})

test_that("wombat file has no header by default", {
  tp  <- tidyped(make_simple_ped())
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)

  pedexport(tp, software = "wombat", file = tmp)

  lines <- readLines(tmp)
  # Wombat reads the pedigree in free format; a header line would be
  # parsed as data, so none is written by default.
  expect_equal(length(lines), nrow(tp))
  expect_false(lines[[1]] == "animal sire dam")
})

test_that("pedexport returns invisible result even when writing file", {
  tp  <- tidyped(make_simple_ped())
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)

  # withVisible to check invisibility
  res <- withVisible(pedexport(tp, software = "blupf90", file = tmp))
  expect_false(res$visible)
  expect_s3_class(res$value, "data.table")
})

# ---------------------------------------------------------------------------
# 5. Failure modes
# ---------------------------------------------------------------------------

test_that("pedexport errors on non-tidyped input without core columns", {
  expect_error(pedexport(data.frame(x = 1:3)), regexp = "tidyped")
})

test_that("pedexport errors on invalid software argument", {
  tp <- tidyped(make_simple_ped())
  expect_error(pedexport(tp, software = "spss"), regexp = "arg")
})

test_that("pedexport errors on invalid file argument", {
  tp <- tidyped(make_simple_ped())
  expect_error(pedexport(tp, file = 123), regexp = "file")
  expect_error(pedexport(tp, file = ""),  regexp = "file")
})

test_that("pedexport errors on invalid header argument", {
  tp <- tidyped(make_simple_ped())
  expect_error(pedexport(tp, header = "yes"), regexp = "header")
})

test_that("blupf90 rejects TAB separators", {
  tp <- tidyped(make_simple_ped())
  expect_error(pedexport(tp, software = "blupf90", sep = "\t"),
               regexp = "TAB")
})

# ---------------------------------------------------------------------------
# 6. Large pedigree smoke test (no error, correct row count)
# ---------------------------------------------------------------------------

test_that("pedexport handles big_family_size_ped without error", {
  tp  <- tidyped(big_family_size_ped)
  out <- pedexport(tp, software = "blupf90")

  expect_equal(nrow(out), nrow(tp))
  expect_true(all(out$IndNum > 0L))
  # Parents always precede offspring
  non_f <- out[SireNum > 0L]
  expect_true(all(non_f$SireNum < non_f$IndNum))
})

# ---------------------------------------------------------------------------
# 7. All-founder pedigree (no parents anywhere)
# ---------------------------------------------------------------------------

test_that("pedexport handles a pedigree of only founders", {
  # tidyped() rejects all-founder input ("All parents are missing"), so the
  # export layer cannot receive one through that path.  Build a complete
  # all-founder tidyped object directly with as_tidyped() instead.
  only_founders <- data.frame(
    Ind     = c("F1", "F2", "F3"),
    Sire    = c(NA,  NA,  NA),
    Dam     = c(NA,  NA,  NA),
    Sex     = c("male", "female", "male"),
    Gen     = c(0L, 0L, 0L),
    IndNum  = c(1L, 2L, 3L),
    SireNum = c(0L, 0L, 0L),
    DamNum  = c(0L, 0L, 0L),
    stringsAsFactors = FALSE
  )
  tp  <- as_tidyped(only_founders)
  out <- pedexport(tp, software = "blupf90")

  expect_equal(nrow(out), 3L)
  expect_true(all(out$SireNum == 0L))
  expect_true(all(out$DamNum  == 0L))
})

# ---------------------------------------------------------------------------
# 8. Optional tidyped columns (addnum = FALSE / addgen = FALSE)
# ---------------------------------------------------------------------------

test_that("numeric exports reconstruct indices when addnum = FALSE", {
  tp  <- tidyped(make_simple_ped(), addnum = FALSE)
  out <- pedexport(tp, software = "blupf90")

  expect_named(out, c("IndNum", "SireNum", "DamNum"))
  expect_type(out$IndNum, "integer")
  expect_equal(out$IndNum, sort(out$IndNum))
  # Parents still precede offspring after reconstruction
  non_f <- out[SireNum > 0L]
  expect_true(all(non_f$SireNum < non_f$IndNum))
})

test_that("asreml export orders founders first when addgen = FALSE", {
  tp  <- tidyped(make_simple_ped(), addgen = FALSE)
  out <- pedexport(tp, software = "asreml")

  founders <- tp[is.na(Sire) & is.na(Dam), Ind]
  first_rows <- head(out$animal, length(founders))
  expect_true(all(founders %in% first_rows))
})

test_that("exports work when both addnum and addgen are FALSE", {
  tp  <- tidyped(make_simple_ped(), addnum = FALSE, addgen = FALSE)
  out <- pedexport(tp, software = "blupf90")

  expect_equal(nrow(out), nrow(tp))
  expect_equal(out$IndNum, seq_len(nrow(out)))
})

test_that("numeric formats reject non-integer missing symbols", {
  tp <- tidyped(make_simple_ped())
  expect_error(pedexport(tp, software = "blupf90", missing = "."),
               regexp = "integer")
  expect_error(pedexport(tp, software = "blupf90", missing = "0"),
               regexp = "integer")

  # Integer-valued numerics are accepted and replace the internal 0 sentinel
  out <- pedexport(tp, software = "blupf90", missing = -99)
  founders <- tp[is.na(Sire) & is.na(Dam), IndNum]
  expect_true(all(out[IndNum %in% founders, SireNum] == -99L))
  expect_true(all(out[IndNum %in% founders, DamNum]  == -99L))
})

# ---------------------------------------------------------------------------
# 9. dmu and sommer formats; numeric ID mapping (xref)
# ---------------------------------------------------------------------------

test_that("pedexport dmu is identical to blupf90", {
  tp  <- tidyped(make_simple_ped())
  b90 <- pedexport(tp, software = "blupf90")
  dmu <- pedexport(tp, software = "dmu")
  expect_equal(b90, dmu)
})

test_that("sommer format returns a character pedigree with NA missing", {
  tp  <- tidyped(make_simple_ped())
  out <- pedexport(tp, software = "sommer")

  expect_s3_class(out, "data.table")
  expect_named(out, c("ID", "Sire", "Dam"))
  expect_type(out$ID, "character")
  # Missing parents are NA (R convention), non-missing keep character IDs
  founders <- tp[is.na(Sire) & is.na(Dam), Ind]
  expect_true(all(out[ID %in% founders, is.na(Sire)]))
  expect_true(all(out[ID %in% founders, is.na(Dam)]))
  expect_equal(out[ID == "C", Sire], "A")
  expect_equal(out[ID == "C", Dam],  "B")
  # Parents precede offspring
  expect_true(which(out$ID == "C") < which(out$ID == "D"))
})

test_that("sommer output pairs with the pedmat relationship matrix", {
  tp  <- tidyped(make_simple_ped())
  out <- pedexport(tp, software = "sommer")
  A   <- pedmat(tp, method = "A")

  # pedmat() rownames are the character IDs in the same order as out$ID,
  # so the matrix can be passed directly to sommer::mmer(..., Gu = A).
  expect_equal(rownames(A), out$ID)
  expect_equal(colnames(A), out$ID)
})

test_that("sommer always codes missing parents as NA regardless of 'missing'", {
  tp  <- tidyped(make_simple_ped())
  out <- pedexport(tp, software = "sommer", missing = "0")

  founders <- tp[is.na(Sire) & is.na(Dam), Ind]
  expect_true(all(out[ID %in% founders, is.na(Sire)]))
})

test_that("numeric exports carry an xref mapping attribute", {
  tp <- tidyped(make_simple_ped())
  for (sw in c("blupf90", "mtdfreml", "dmu", "numeric")) {
    out <- pedexport(tp, software = sw)
    xr  <- attr(out, "xref")
    expect_s3_class(xr, "data.table")
    expect_named(xr, c("IndNum", "Ind"))
    expect_equal(xr$IndNum, seq_len(nrow(tp)))
    expect_equal(xr$Ind, tp$Ind)
  }
})

test_that("character formats carry no xref attribute", {
  tp <- tidyped(make_simple_ped())
  expect_null(attr(pedexport(tp, software = "asreml"), "xref"))
  expect_null(attr(pedexport(tp, software = "wombat"), "xref"))
  expect_null(attr(pedexport(tp, software = "sommer"), "xref"))
})

test_that("writing a numeric file also writes the .xref mapping file", {
  tp  <- tidyped(make_simple_ped())
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(c(tmp, paste0(tmp, ".xref"))), add = TRUE)

  pedexport(tp, software = "blupf90", file = tmp)

  xr_file <- paste0(tmp, ".xref")
  expect_true(file.exists(xr_file))
  xr <- read.table(xr_file, header = FALSE, stringsAsFactors = FALSE)
  expect_equal(ncol(xr), 2L)
  # Row 1 maps numeric 1 to the first individual
  expect_equal(xr[1, 1], 1L)
  expect_equal(xr[1, 2], tp$Ind[1])
})
