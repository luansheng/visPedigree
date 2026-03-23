library(testthat)
library(visPedigree)
library(data.table)

# ===========================================================================
# Test Suite for pediv()
# ===========================================================================

# ---- Fixtures (reuse theory pedigree from test-pedcontrib.R) ----
make_theory_ped <- function() {
  ped_dt <- data.table(
    Ind  = c("F1", "F2", "F3", "F4", "S1", "D1", "R1", "R2", "R3", "R4"),
    Sire = c(NA,   NA,   NA,   NA,  "F1", "F3", "S1", "S1", "S1", "S1"),
    Dam  = c(NA,   NA,   NA,   NA,  "F2", "F4", "D1", "D1", "D1", "F4")
  )
  suppressMessages(tidyped(ped_dt))
}

make_small_ped <- function() {
  suppressMessages(tidyped(small_ped))
}


# ==========================================================================
# Group 1: Return structure
# ==========================================================================

test_that("pediv returns correct class and list structure", {
  tp <- make_theory_ped()
  div <- suppressMessages(pediv(tp, reference = c("R1", "R2", "R3", "R4")))

  expect_s3_class(div, "pediv")
  expect_true(is.list(div))
  expect_true(all(c("summary", "founders", "ancestors") %in% names(div)))
})

test_that("pediv summary is a single-row data.table with expected columns", {
  tp <- make_theory_ped()
  div <- suppressMessages(pediv(tp, reference = c("R1", "R2", "R3", "R4")))

  expected_cols <- c("NRef", "NFounder", "feH", "fe",
                     "NAncestor", "faH", "fa", "fafe",
                     "NeCoancestry", "NeInbreeding", "NeDemographic")
  expect_true(all(expected_cols %in% names(div$summary)))
  expect_equal(nrow(div$summary), 1L)
})


# ==========================================================================
# Group 2: f_e and f_a values match pedcontrib directly
# ==========================================================================

test_that("pediv f_e matches pedcontrib for theory pedigree", {
  tp <- make_theory_ped()
  ref <- c("R1", "R2", "R3", "R4")

  div  <- suppressMessages(pediv(tp, reference = ref))
  cont <- suppressMessages(pedcontrib(tp, reference = ref, mode = "both"))

  expect_equal(div$summary$fe, cont$summary$f_e, tolerance = 1e-10)
})

test_that("pediv f_a matches pedcontrib for theory pedigree", {
  tp <- make_theory_ped()
  ref <- c("R1", "R2", "R3", "R4")

  div  <- suppressMessages(pediv(tp, reference = ref))
  cont <- suppressMessages(pedcontrib(tp, reference = ref, mode = "both"))

  expect_equal(div$summary$fa, cont$summary$f_a, tolerance = 1e-10)
})

test_that("pediv f_e and f_a match known values for theory pedigree", {
  tp  <- make_theory_ped()
  div <- suppressMessages(pediv(tp, reference = c("R1", "R2", "R3", "R4")))

  expect_equal(div$summary$fe, 1 / (0.25^2 + 0.25^2 + 0.1875^2 + 0.3125^2),
               tolerance = 1e-6)
  expect_equal(div$summary$fa, 1 / (0.5^2 + 0.375^2 + 0.125^2),
               tolerance = 1e-6)
})

test_that("pediv small_ped: f_e and f_a match known values", {
  tp  <- make_small_ped()
  div <- suppressMessages(pediv(tp, reference = c("Z1", "Z2", "Y", "X")))

  expect_equal(div$summary$fe, 6.585209, tolerance = 1e-4)
  expect_equal(div$summary$fa, 2.666667, tolerance = 1e-4)
})


# ==========================================================================
# Group 3: f_a_fe_ratio
# ==========================================================================

test_that("fafe equals fa / fe", {
  tp  <- make_theory_ped()
  div <- suppressMessages(pediv(tp, reference = c("R1", "R2", "R3", "R4")))

  expect_equal(div$summary$fafe,
               div$summary$fa / div$summary$fe,
               tolerance = 1e-6)
})

test_that("fafe is between 0 and 1", {
  tp  <- make_small_ped()
  div <- suppressMessages(pediv(tp, reference = c("Z1", "Z2", "Y", "X")))

  expect_true(div$summary$fafe > 0)
  expect_true(div$summary$fafe <= 1)
})


# ==========================================================================
# Group 4: Ne values are numeric (or NA, never an error)
# ==========================================================================

test_that("Ne columns are numeric (not character or list)", {
  tp  <- make_theory_ped()
  div <- suppressMessages(pediv(tp, reference = c("R1", "R2", "R3", "R4")))

  expect_true(is.numeric(div$summary$NeCoancestry))
  expect_true(is.numeric(div$summary$NeInbreeding))
  expect_true(is.numeric(div$summary$NeDemographic))
})

test_that("Ne_coancestry is positive when computable", {
  tp  <- make_small_ped()
  div <- suppressMessages(pediv(tp, reference = c("Z1", "Z2", "Y", "X")))

  # coancestry method should produce a positive Ne for a real pedigree
  if (!is.na(div$summary$NeCoancestry)) {
    expect_true(div$summary$NeCoancestry > 0)
  }
})


# ==========================================================================
# Group 5: n_ref propagates correctly
# ==========================================================================

test_that("n_ref reflects length of reference argument", {
  tp  <- make_theory_ped()
  ref <- c("R1", "R2", "R3", "R4")
  div <- suppressMessages(pediv(tp, reference = ref))

  expect_equal(div$summary$NRef, length(ref))
})


# ==========================================================================
# Group 6: founders / ancestors tables are propagated from pedcontrib
# ==========================================================================

test_that("founders table matches pedcontrib", {
  tp  <- make_theory_ped()
  ref <- c("R1", "R2", "R3", "R4")

  div  <- suppressMessages(pediv(tp, reference = ref, top = 100))
  cont <- suppressMessages(pedcontrib(tp, reference = ref, mode = "both", top = 100))

  expect_equal(div$founders, cont$founders)
})

test_that("ancestors table matches pedcontrib", {
  tp  <- make_theory_ped()
  ref <- c("R1", "R2", "R3", "R4")

  div  <- suppressMessages(pediv(tp, reference = ref, top = 100))
  cont <- suppressMessages(pedcontrib(tp, reference = ref, mode = "both", top = 100))

  expect_equal(div$ancestors, cont$ancestors)
})


# ==========================================================================
# Group 7: Input validation
# ==========================================================================

test_that("pediv rejects non-tidyped input", {
  expect_error(pediv(data.frame(a = 1)), "tidyped")
})

test_that("pediv feH matches pedcontrib f_e_H", {
  tp  <- make_theory_ped()
  ref <- c("R1", "R2", "R3", "R4")
  div  <- suppressMessages(pediv(tp, reference = ref))
  cont <- suppressMessages(pedcontrib(tp, reference = ref, mode = "both"))

  expect_equal(div$summary$feH, cont$summary$f_e_H, tolerance = 1e-10)
})

test_that("pediv faH matches pedcontrib f_a_H", {
  tp  <- make_theory_ped()
  ref <- c("R1", "R2", "R3", "R4")
  div  <- suppressMessages(pediv(tp, reference = ref))
  cont <- suppressMessages(pedcontrib(tp, reference = ref, mode = "both"))

  expect_equal(div$summary$faH, cont$summary$f_a_H, tolerance = 1e-10)
})

test_that("print.pediv runs without error", {
  tp  <- make_theory_ped()
  div <- suppressMessages(pediv(tp, reference = c("R1", "R2", "R3", "R4")))
  expect_output(print(div), "Genetic Diversity Summary")
  expect_output(print(div), "fe\\(H\\)")
  expect_output(print(div), "fa\\(H\\)")
  expect_output(print(div), "Ne")
})


# ==========================================================================
# Group 8: GeneDiv
# ==========================================================================

test_that("GeneDiv column is present in summary", {
  tp  <- make_theory_ped()
  div <- suppressMessages(pediv(tp, reference = c("R1", "R2", "R3", "R4")))
  expect_true("GeneDiv" %in% names(div$summary))
})

test_that("GeneDiv equals 1 - MeanCoan", {
  tp  <- make_small_ped()
  div <- suppressMessages(pediv(tp, reference = c("Z1", "Z2", "Y", "X")))
  if (!is.na(div$summary$MeanCoan) && !is.na(div$summary$GeneDiv)) {
    expect_equal(div$summary$GeneDiv, 1 - div$summary$MeanCoan, tolerance = 1e-12)
  }
})

test_that("GeneDiv is in [0, 1]", {
  tp  <- make_small_ped()
  div <- suppressMessages(pediv(tp, reference = c("Z1", "Z2", "Y", "X")))
  if (!is.na(div$summary$GeneDiv)) {
    expect_true(div$summary$GeneDiv >= 0)
    expect_true(div$summary$GeneDiv <= 1)
  }
})

test_that("print.pediv displays GeneDiv", {
  tp  <- make_small_ped()
  div <- suppressMessages(pediv(tp, reference = c("Z1", "Z2", "Y", "X")))
  expect_output(print(div), "GeneDiv")
})
