library(testthat)
library(visPedigree)
library(data.table)

make_pedrel_test_pedigree <- function() {
  tidyped(data.table(
    Ind = c("A", "B", "C", "D", "E", "F", "G", "H"),
    Sire = c(NA, NA, "A", "A", "C", "C", "E", "E"),
    Dam = c(NA, NA, "B", "B", "D", "D", "F", "F"),
    Cohort = rep(1:4, each = 2)
  ))
}

explicit_group_summary <- function(group, ped, scale, reference = NULL) {
  ids <- ped[Cohort == group, Ind]
  if (!is.null(reference)) ids <- ids[ids %in% reference]
  A <- pedmat(ped, method = "A", sparse = FALSE)
  idx <- match(ids, rownames(A))
  A_group <- A[idx, idx, drop = FALSE]

  if (scale == "relationship") {
    (sum(A_group) - sum(diag(A_group))) /
      (length(idx) * (length(idx) - 1))
  } else {
    sum(A_group) / (2 * length(idx)^2)
  }
}

test_that("pedrel matrix-free summaries match an explicit A matrix", {
  ped <- make_pedrel_test_pedigree()

  rel <- pedrel(ped, by = "Cohort", scale = "relationship")
  coan <- pedrel(ped, by = "Cohort", scale = "coancestry")

  expected_rel <- vapply(
    rel$Cohort, explicit_group_summary, numeric(1),
    ped = ped, scale = "relationship"
  )
  expected_coan <- vapply(
    coan$Cohort, explicit_group_summary, numeric(1),
    ped = ped, scale = "coancestry"
  )

  expect_equal(rel$MeanRel, expected_rel, tolerance = 1e-12)
  expect_equal(coan$MeanCoan, expected_coan, tolerance = 1e-12)
  expect_true(all(rel$Status == "ok"))
  expect_true(all(coan$Status == "ok"))
})

test_that("pedrel reference summaries match the selected explicit submatrix", {
  ped <- make_pedrel_test_pedigree()
  reference <- c("C", "D", "G", "H")

  rel <- suppressWarnings(
    pedrel(ped, by = "Cohort", reference = reference)
  )
  coan <- suppressWarnings(
    pedrel(ped, by = "Cohort", reference = reference, scale = "coancestry")
  )

  valid_rel <- rel[Status == "ok"]
  valid_coan <- coan[Status == "ok"]
  expected_rel <- vapply(
    valid_rel$Cohort, explicit_group_summary, numeric(1),
    ped = ped, scale = "relationship", reference = reference
  )
  expected_coan <- vapply(
    valid_coan$Cohort, explicit_group_summary, numeric(1),
    ped = ped, scale = "coancestry", reference = reference
  )

  expect_equal(valid_rel$MeanRel, expected_rel, tolerance = 1e-12)
  expect_equal(valid_coan$MeanCoan, expected_coan, tolerance = 1e-12)
})

test_that("pedrel handles a reference consisting only of founders", {
  ped <- make_pedrel_test_pedigree()
  reference <- c("A", "B")

  rel <- suppressWarnings(
    pedrel(ped, by = "Cohort", reference = reference)
  )
  coan <- suppressWarnings(
    pedrel(ped, by = "Cohort", reference = reference, scale = "coancestry")
  )

  expect_equal(rel[Cohort == 1, MeanRel], 0)
  expect_equal(coan[Cohort == 1, MeanCoan], 0.25)
  expect_equal(rel[Cohort == 1, Status], "ok")
  expect_equal(coan[Cohort == 1, Status], "ok")
})

test_that("pedrel matches explicit A for one-parent and selfing pedigrees", {
  one_parent <- tidyped(data.table(
    Ind = c("A", "B", "C"),
    Sire = c(NA, "A", "A"),
    Dam = c(NA, NA, NA),
    Cohort = c(NA, 1L, 1L)
  ))
  selfed <- suppressMessages(tidyped(
    data.table(
      Ind = c("A", "B", "C", "D", "E"),
      Sire = c(NA, NA, "A", "C", "C"),
      Dam = c(NA, NA, "B", "C", "C"),
      Cohort = c(NA, NA, NA, 1L, 1L)
    ),
    selfing = TRUE
  ))

  for (ped in list(one_parent, selfed)) {
    rel <- pedrel(ped, by = "Cohort")
    coan <- pedrel(ped, by = "Cohort", scale = "coancestry")

    expect_equal(
      rel$MeanRel,
      explicit_group_summary(1L, ped, "relationship"),
      tolerance = 1e-12
    )
    expect_equal(
      coan$MeanCoan,
      explicit_group_summary(1L, ped, "coancestry"),
      tolerance = 1e-12
    )
  }
})

test_that("pedrel no longer constructs or scans a relationship matrix", {
  ped <- make_pedrel_test_pedigree()
  local_mocked_bindings(
    pedmat = function(...) stop("pedmat() must not be called"),
    .package = "visPedigree"
  )

  rel <- pedrel(ped, by = "Cohort")
  coan <- pedrel(ped, by = "Cohort", scale = "coancestry", compact = TRUE)

  expect_true(all(rel$Status == "ok"))
  expect_true(all(coan$Status == "ok"))
})

test_that("pedrel processes more groups than one product batch", {
  n_groups <- 40L
  raw <- data.table(
    Ind = c("Sire", "Dam", sprintf("F%03d", seq_len(2L * n_groups))),
    Sire = c(NA, NA, rep("Sire", 2L * n_groups)),
    Dam = c(NA, NA, rep("Dam", 2L * n_groups)),
    Cohort = c(NA_integer_, NA_integer_, rep(seq_len(n_groups), each = 2L))
  )
  ped <- tidyped(raw)

  rel <- pedrel(ped, by = "Cohort")
  coan <- pedrel(ped, by = "Cohort", scale = "coancestry")

  expect_equal(rel$MeanRel, rep(0.5, n_groups))
  expect_equal(coan$MeanCoan, rep(0.375, n_groups))
  expect_true(all(rel$Status == "ok"))
  expect_true(all(coan$Status == "ok"))
})

test_that("pedrel validates compact", {
  ped <- make_pedrel_test_pedigree()

  expect_error(pedrel(ped, by = "Cohort", compact = NA), "compact")
})
