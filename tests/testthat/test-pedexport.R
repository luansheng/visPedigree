# Tests for pedexport()

library(testthat)
library(data.table)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_simple_ped <- function() {
  # A -> C, B -> C, C -> D (simple 3-generation chain with one inbred offspring)
  data.frame(
    Ind  = c("A", "B", "C", "D"),
    Sire = c(NA,  NA,  "A", "C"),
    Dam  = c(NA,  NA,  "B", "C"),
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

test_that("pedexport wombat is identical to blupf90", {
  tp   <- tidyped(make_simple_ped())
  b90  <- pedexport(tp, software = "blupf90")
  womb <- pedexport(tp, software = "wombat")
  expect_equal(b90, womb)
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

test_that("pedexport writes a file with header for asreml", {
  tp  <- tidyped(make_simple_ped())
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)

  pedexport(tp, software = "asreml", file = tmp)

  lines <- readLines(tmp)
  # First line is header
  expect_equal(lines[[1]], "animal sire dam")
  # Total lines = header + n individuals
  expect_equal(length(lines), nrow(tp) + 1L)
})

test_that("pedexport respects header = FALSE for asreml", {
  tp  <- tidyped(make_simple_ped())
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)

  pedexport(tp, software = "asreml", file = tmp, header = FALSE)

  lines <- readLines(tmp)
  expect_equal(length(lines), nrow(tp))
  # First line must NOT be the header
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
  only_founders <- data.frame(
    Ind  = c("F1", "F2", "F3"),
    Sire = c(NA, NA, NA),
    Dam  = c(NA, NA, NA),
    stringsAsFactors = FALSE
  )
  tp  <- tidyped(only_founders)
  out <- pedexport(tp, software = "blupf90")

  expect_equal(nrow(out), 3L)
  expect_true(all(out$SireNum == 0L))
  expect_true(all(out$DamNum  == 0L))
})
