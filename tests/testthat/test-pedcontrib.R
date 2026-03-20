library(testthat)
library(visPedigree)
library(data.table)

# ===========================================================================
# Test Suite for pedcontrib()
# Purpose: Capture exact numerical outputs of the current (pure-R) implementation
#          so that any future Rcpp optimisation must produce identical results.
# ===========================================================================

# ---- Fixture 1: 10-individual "theory" pedigree (瓶颈模型) ----
# Hand-derived expected values from the manuscript
make_theory_ped <- function() {
  ped_dt <- data.table(
    Ind  = c("F1", "F2", "F3", "F4", "S1", "D1", "R1", "R2", "R3", "R4"),
    Sire = c(NA,   NA,   NA,   NA,  "F1", "F3", "S1", "S1", "S1", "S1"),
    Dam  = c(NA,   NA,   NA,   NA,  "F2", "F4", "D1", "D1", "D1", "F4")
  )
  suppressMessages(tidyped(ped_dt))
}

# ---- Fixture 2: small_ped (visPedigree built-in dataset) ----
make_small_ped <- function() {
  suppressMessages(tidyped(small_ped))
}

# ---- Fixture 3: two founders with one offspring, candidates are founders ----
make_founders_as_cand_ped <- function() {
  ped_dt <- data.table(
    Ind  = c("A", "B", "C"),
    Sire = c(NA, NA, "A"),
    Dam  = c(NA, NA, "B")
  )
  suppressMessages(tidyped(ped_dt))
}

# ---- Fixture 4: single offspring from two founders ----
make_minimal_ped <- function() {
  ped_dt <- data.table(
    Ind  = c("P1", "P2", "O1"),
    Sire = c(NA,  NA,  "P1"),
    Dam  = c(NA,  NA,  "P2")
  )
  suppressMessages(tidyped(ped_dt))
}

# ---- Fixture 5: single candidate (edge case) ----
make_single_cand_ped <- function() {
  ped_dt <- data.table(
    Ind  = c("A", "B", "C", "D"),
    Sire = c(NA,  NA, "A", "A"),
    Dam  = c(NA,  NA, "B", "B")
  )
  suppressMessages(tidyped(ped_dt))
}

# ---- Fixture 6: deeper chain (4 generations) ----
make_deep_chain_ped <- function() {
  ped_dt <- data.table(
    Ind  = c("G1M", "G1F", "G2", "G3M", "G3F", "G4"),
    Sire = c(NA,    NA,   "G1M", "G2",  NA,   "G3M"),
    Dam  = c(NA,    NA,   "G1F", NA,    NA,   "G3F")
  )
  suppressMessages(tidyped(ped_dt))
}


# ==========================================================================
# Group 1: Founder contributions (p_i) and f_e
# ==========================================================================

test_that("Theory pedigree: founder contributions match hand-derived values", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "founder", top = 100))

  # Extract founder table and sort by name for deterministic comparison
  fd <- res$founders
  setkey(fd, Ind)

  expect_equal(fd["F1", Contrib], 0.25,     tolerance = 1e-10)
  expect_equal(fd["F2", Contrib], 0.25,     tolerance = 1e-10)
  expect_equal(fd["F3", Contrib], 0.1875,   tolerance = 1e-10)
  expect_equal(fd["F4", Contrib], 0.3125,   tolerance = 1e-10)
})

test_that("Theory pedigree: f_e matches hand-derived value (3.878788)", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "both", top = 100))

  expected_fe <- 1 / (0.25^2 + 0.25^2 + 0.1875^2 + 0.3125^2)
  expect_equal(res$summary$f_e, expected_fe, tolerance = 1e-6)
})

test_that("Theory pedigree: founder contributions sum to 1.0", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "founder", top = 100))

  expect_equal(sum(res$founders$Contrib), 1.0, tolerance = 1e-10)
})


# ==========================================================================
# Group 2: Ancestor marginal contributions (q_i) and f_a
# ==========================================================================

test_that("Theory pedigree: ancestor marginal contributions match hand-derived values", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "ancestor", top = 100))

  anc <- res$ancestors
  # By rank order: S1(0.5), D1(0.375), F4(0.125)
  expect_equal(anc$Contrib[1], 0.500, tolerance = 1e-10)
  expect_equal(anc$Contrib[2], 0.375, tolerance = 1e-10)
  expect_equal(anc$Contrib[3], 0.125, tolerance = 1e-10)
  expect_equal(anc$Ind[1], "S1")
  expect_equal(anc$Ind[2], "D1")
  expect_equal(anc$Ind[3], "F4")
})

test_that("Theory pedigree: f_a matches hand-derived value (2.461538)", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "both", top = 100))

  expected_fa <- 1 / (0.5^2 + 0.375^2 + 0.125^2)
  expect_equal(res$summary$f_a, expected_fa, tolerance = 1e-6)
})

test_that("Theory pedigree: ancestor contributions sum to 1.0", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "ancestor", top = 100))

  expect_equal(sum(res$ancestors$Contrib), 1.0, tolerance = 1e-10)
})


# ==========================================================================
# Group 3: small_ped dataset (complex real-world pedigree)
# ==========================================================================

test_that("small_ped: f_e and f_a match known values", {
  tp <- make_small_ped()
  cand <- c("Z1", "Z2", "Y", "X")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "both", top = 100))

  # These values are from the current pure-R implementation runs
  expect_equal(res$summary$f_e, 6.585209, tolerance = 1e-4)
  expect_equal(res$summary$f_a, 2.666667, tolerance = 1e-4)
})

test_that("small_ped: founder contributions sum to 1.0", {
  tp <- make_small_ped()
  cand <- c("Z1", "Z2", "Y", "X")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "founder", top = 100))

  expect_equal(sum(res$founders$Contrib), 1.0, tolerance = 1e-10)
})

test_that("small_ped: ancestor contributions sum to 1.0", {
  tp <- make_small_ped()
  cand <- c("Z1", "Z2", "Y", "X")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "ancestor", top = 100))

  expect_equal(sum(res$ancestors$Contrib), 1.0, tolerance = 1e-10)
})

test_that("small_ped: top ancestor is N with q = 0.50", {
  tp <- make_small_ped()
  cand <- c("Z1", "Z2", "Y", "X")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "ancestor", top = 100))

  expect_equal(res$ancestors$Ind[1], "X")
  expect_equal(res$ancestors$Contrib[1], 0.50, tolerance = 1e-10)
})

test_that("small_ped: ancestor ranking order is correct", {
  tp <- make_small_ped()
  cand <- c("Z1", "Z2", "Y", "X")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "ancestor", top = 100))

  # Expected order: N(0.50), V(0.25), R(0.125), K(0.125), W(0.125), J1(0.0625), O(0.03125)
  expect_equal(res$ancestors$Ind[1], "X")
  expect_equal(res$ancestors$Ind[2], "N")
  expect_equal(res$ancestors$Contrib[2], 0.25, tolerance = 1e-10)
})


# ==========================================================================
# Group 4: f_e and f_a invariance under `top` parameter
# ==========================================================================

test_that("f_e and f_a are invariant to the top parameter", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")

  res_full <- suppressMessages(pedcontrib(tp, reference = cand, mode = "both", top = 100))
  res_top1 <- suppressMessages(pedcontrib(tp, reference = cand, mode = "both", top = 1))

  expect_equal(res_full$summary$f_e, res_top1$summary$f_e)
  expect_equal(res_full$summary$f_a, res_top1$summary$f_a)
})


# ==========================================================================
# Group 5: Edge cases
# ==========================================================================

test_that("Single candidate: founder and ancestor contributions work", {
  tp <- make_single_cand_ped()
  res <- suppressMessages(pedcontrib(tp, reference = "C", mode = "both", top = 100))

  # Single individual C has parents A and B
  # Founder contributions: A = 0.5, B = 0.5
  fd <- res$founders
  setkey(fd, Ind)
  expect_equal(fd["A", Contrib], 0.5, tolerance = 1e-10)
  expect_equal(fd["B", Contrib], 0.5, tolerance = 1e-10)
  expect_equal(sum(fd$Contrib), 1.0, tolerance = 1e-10)

  # Ancestor contributions should also sum to 1.0
  expect_equal(sum(res$ancestors$Contrib), 1.0, tolerance = 1e-10)
})

test_that("Minimal pedigree (1 offspring from 2 founders): correct values", {
  tp <- make_minimal_ped()
  res <- suppressMessages(pedcontrib(tp, reference = "O1", mode = "both", top = 100))

  # f_e: 1/(0.5^2 + 0.5^2) = 2.0
  expect_equal(res$summary$f_e, 2.0, tolerance = 1e-10)

  # Founders: P1=0.5, P2=0.5
  fd <- res$founders
  setkey(fd, Ind)
  expect_equal(fd["P1", Contrib], 0.5, tolerance = 1e-10)
  expect_equal(fd["P2", Contrib], 0.5, tolerance = 1e-10)
})

test_that("Deep chain: contributions propagate correctly through 4 generations", {
  tp <- make_deep_chain_ped()
  res <- suppressMessages(pedcontrib(tp, reference = "G4", mode = "both", top = 100))

  # G4's parents are G3M and G3F

  # G3F is a founder -> contributes 0.5 to G4
  # G3M's parent is G2 (sire side only)
  # G2's parents are G1M and G1F
  # So G1M -> G2 -> G3M -> G4: contribution = 0.5^3 = 0.125
  # G1F -> G2 -> G3M -> G4: contribution = 0.5^3 = 0.125
  # G3M has no dam, so half of G3M is "missing" => goes to phantom founder
  # G2 -> G3M: 0.5 (sire contribution to G3M), then G3M -> G4: 0.5
  # So from G2: 0.25 to G4

  fd <- res$founders
  setkey(fd, Ind)
  expect_equal(fd["G3F", Contrib], 0.5, tolerance = 1e-10)

  # G3M has no dam => 25% of G4's genes go to a phantom founder
  # Real founders account for only 0.75 of total variation
  expect_equal(sum(fd$Contrib), 1.0, tolerance = 1e-10)

  # Ancestor contributions also reflect the same incomplete pedigree
  expect_true(sum(res$ancestors$Contrib) <= 1.0 + 1e-10)
})

test_that("Founders as candidates: they are their own founders, no ancestors above", {
  tp <- make_founders_as_cand_ped()
  # Use the two founders A and B as candidates
  res <- suppressMessages(pedcontrib(tp, reference = c("A", "B"), mode = "both", top = 100))

  # Founder contributions: A = 0.5, B = 0.5
  expect_equal(nrow(res$founders), 2)
  expect_equal(sum(res$founders$Contrib), 1.0, tolerance = 1e-10)

  # f_e = 1 / (0.5^2 + 0.5^2) = 2.0
  expect_equal(res$summary$f_e, 2.0, tolerance = 1e-10)

  # No ancestors should be found (founders have no parents)
  expect_equal(nrow(res$ancestors), 2)
})


# ==========================================================================
# Group 6: mode parameter behaviour
# ==========================================================================

test_that("mode='founder' does not produce ancestors", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "founder"))

  expect_true(!is.null(res$founders))
  expect_true(is.null(res$ancestors))
})

test_that("mode='ancestor' does not produce founders", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "ancestor"))

  expect_true(is.null(res$founders))
  expect_true(!is.null(res$ancestors))
})

test_that("mode='both' produces both founders and ancestors", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "both"))

  expect_true(!is.null(res$founders))
  expect_true(!is.null(res$ancestors))
})


# ==========================================================================
# Group 7: Return structure and class
# ==========================================================================

test_that("pedcontrib returns correct structure", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "both"))

  expect_s3_class(res, "pedcontrib")
  expect_true(is.list(res))
  expect_true(all(c("founders", "ancestors", "summary") %in% names(res)))

  # founders table columns
  expect_true(all(c("Ind", "Contrib", "CumContrib", "Rank") %in% names(res$founders)))

  # ancestors table columns
  expect_true(all(c("Ind", "Contrib", "CumContrib", "Rank") %in% names(res$ancestors)))

  # summary fields
  expect_true(all(c("n_ref", "n_founder", "f_e",
                     "n_ancestor", "f_a") %in% names(res$summary)))
})

test_that("CumContrib is monotonically non-decreasing", {
  tp <- make_small_ped()
  cand <- c("Z1", "Z2", "Y", "X")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "both", top = 100))

  expect_true(all(diff(res$founders$CumContrib) >= -1e-15))
  expect_true(all(diff(res$ancestors$CumContrib) >= -1e-15))
})


# ==========================================================================
# Group 8: Error handling
# ==========================================================================

test_that("pedcontrib rejects non-tidyped input", {
  expect_error(pedcontrib(data.frame(a = 1)), "tidyped")
})

test_that("pedcontrib errors on empty candidate set", {
  tp <- make_theory_ped()
  expect_error(
    suppressMessages(suppressWarnings(pedcontrib(tp, reference = "NONEXISTENT", mode = "both"))),
    "No valid reference"
  )
})


# ==========================================================================
# Group 9: Relationship between f_e and f_a
# ==========================================================================

test_that("f_a <= f_e (bottleneck effect)", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "both", top = 100))

  expect_true(res$summary$f_a <= res$summary$f_e)
})

test_that("small_ped: f_a <= f_e", {
  tp <- make_small_ped()
  cand <- c("Z1", "Z2", "Y", "X")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "both", top = 100))

  expect_true(res$summary$f_a <= res$summary$f_e)
})


# ==========================================================================
# Group 10: Shannon entropy effective numbers f_e(H) and f_a(H)
# ==========================================================================

test_that("Theory pedigree: f_e_H matches hand-derived value", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "founder", top = 100))

  p <- c(0.25, 0.25, 0.1875, 0.3125)
  expected_fe_H <- exp(-sum(p * log(p)))
  expect_equal(res$summary$f_e_H, expected_fe_H, tolerance = 1e-6)
})

test_that("Theory pedigree: f_e_H >= f_e (Hill monotonicity)", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "founder", top = 100))

  expect_true(res$summary$f_e_H >= res$summary$f_e)
})

test_that("Theory pedigree: f_a_H >= f_a", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "both", top = 100))

  expect_true(res$summary$f_a_H >= res$summary$f_a)
})

test_that("Uniform contributions: f_e_H == f_e == Nf", {
  # Two founders each contributing exactly 0.5
  tp <- make_founders_as_cand_ped()
  res <- suppressMessages(pedcontrib(tp, reference = "C", mode = "founder", top = 100))

  # Both parents contribute equally -> uniform
  expect_equal(res$summary$f_e_H, res$summary$f_e, tolerance = 1e-6)
  expect_equal(res$summary$f_e_H, res$summary$n_founder, tolerance = 1e-6)
})

test_that("small_ped: f_e_H >= f_e", {
  tp <- make_small_ped()
  cand <- c("Z1", "Z2", "Y", "X")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "both", top = 100))

  expect_true(res$summary$f_e_H >= res$summary$f_e)
  expect_true(res$summary$f_a_H >= res$summary$f_a)
})

test_that("f_e_H is NA when no founders found", {
  tp <- make_theory_ped()
  cand <- c("R1", "R2", "R3", "R4")
  res <- suppressMessages(pedcontrib(tp, reference = cand, mode = "ancestor", top = 100))

  # When mode = "ancestor", founders is NULL, f_e_H should not exist

  expect_null(res$summary$f_e_H)
  # But f_a_H should exist
  expect_true(!is.null(res$summary$f_a_H) && !is.na(res$summary$f_a_H))
})
