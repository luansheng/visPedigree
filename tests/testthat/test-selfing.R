# Tests for selfing support (Issue #10: plant breeding monoecious pedigrees)

test_that("selfing parameter validation", {
  skip_if_not_installed("visPedigree")
  
  ped <- data.frame(
    Ind = c("A", "B", "C"),
    Sire = c(NA, NA, "A"),
    Dam = c(NA, NA, "B")
  )
  
  # selfing must be logical

  expect_error(tidyped(ped, selfing = 1), "only is assigned using TRUE or FALSE")
  expect_error(tidyped(ped, selfing = "yes"), "only is assigned using TRUE or FALSE")
  
  # selfing = FALSE is default and works
  res <- tidyped(ped, selfing = FALSE)
  expect_s3_class(res, "tidyped")
})

test_that("selfing = TRUE allows same individual as both Sire and Dam", {
  skip_if_not_installed("visPedigree")
  
  # Simple selfing: A selfs to produce C
  ped <- data.frame(
    Ind = c("A", "B", "C"),
    Sire = c(NA, NA, "A"),
    Dam = c(NA, NA, "A"),  # A is both Sire and Dam
    stringsAsFactors = FALSE
  )
  
  # Should error without selfing
  expect_error(tidyped(ped), "Sex conflict detected")
  
  # Should succeed with selfing = TRUE
  expect_message(tidyped(ped, selfing = TRUE), "Selfing mode")
  res <- suppressMessages(tidyped(ped, selfing = TRUE))
  
  expect_s3_class(res, "tidyped")
  expect_equal(res[Ind == "A", Sex], "monoecious")
  expect_true(isTRUE(attr(res, "selfing")))
  expect_equal(attr(res, "bisexual_parents"), "A")
})

test_that("selfing with multiple monoecious individuals", {
  skip_if_not_installed("visPedigree")
  
  # Plant pedigree: P1 and P2 are both monoecious
  ped <- data.frame(
    Ind = c("P1", "P2", "C1", "C2", "C3"),
    Sire = c(NA, NA, "P1", "P2", "P1"),
    Dam = c(NA, NA, "P2", "P1", "P1"),  # P1 and P2 both appear as Sire and Dam
    stringsAsFactors = FALSE
  )
  
  res <- suppressMessages(tidyped(ped, selfing = TRUE))
  
  expect_equal(res[Ind == "P1", Sex], "monoecious")
  expect_equal(res[Ind == "P2", Sex], "monoecious")
  expect_equal(sort(attr(res, "bisexual_parents")), c("P1", "P2"))
})

test_that("selfing: non-monoecious individuals get normal sex inference", {
  skip_if_not_installed("visPedigree")
  
  # Mixed pedigree: P1 is monoecious, M1 is only a Sire, F1 is only a Dam
  ped <- data.frame(
    Ind = c("P1", "M1", "F1", "C1", "C2", "C3"),
    Sire = c(NA, NA, NA, "P1", "M1", "P1"),
    Dam = c(NA, NA, NA, "P1", "F1", "F1"),
    stringsAsFactors = FALSE
  )
  
  res <- suppressMessages(tidyped(ped, selfing = TRUE))
  
  expect_equal(res[Ind == "P1", Sex], "monoecious")
  expect_equal(res[Ind == "M1", Sex], "male")
  expect_equal(res[Ind == "F1", Sex], "female")
})

test_that("selfing: explicit sex annotation conflict for monoecious", {
  skip_if_not_installed("visPedigree")
  
  # A is both Sire and Dam but explicitly annotated as "male" -> should error
  ped <- data.frame(
    Ind = c("A", "B", "C"),
    Sire = c(NA, NA, "A"),
    Dam = c(NA, NA, "A"),
    Sex = c("male", NA, NA),
    stringsAsFactors = FALSE
  )
  
  expect_error(
    tidyped(ped, selfing = TRUE),
    "Sex annotation conflicts for monoecious"
  )
})

test_that("selfing: explicit monoecious annotation accepted", {
  skip_if_not_installed("visPedigree")
  
  # A is both Sire and Dam, explicitly annotated as "monoecious" -> should work
  ped <- data.frame(
    Ind = c("A", "B", "C"),
    Sire = c(NA, NA, "A"),
    Dam = c(NA, NA, "A"),
    Sex = c("monoecious", NA, NA),
    stringsAsFactors = FALSE
  )
  
  res <- suppressMessages(tidyped(ped, selfing = TRUE))
  expect_equal(res[Ind == "A", Sex], "monoecious")
})

test_that("selfing: summary statistics count monoecious correctly", {
  skip_if_not_installed("visPedigree")
  
  ped <- data.frame(
    Ind = c("P1", "P2", "M1", "C1", "C2"),
    Sire = c(NA, NA, NA, "P1", "M1"),
    Dam = c(NA, NA, NA, "P1", "P2"),
    stringsAsFactors = FALSE
  )
  
  res <- suppressMessages(tidyped(ped, selfing = TRUE))
  s <- summary(res)
  
  expect_equal(s$n_monoecious, 1)  # P1
  expect_equal(s$n_male, 1)        # M1
  expect_equal(s$n_female, 1)      # P2 (only Dam)
  
  # n_parents should deduplicate P1 (which is both sire and dam)
  # Sires: P1, M1 -> 2 unique sires
  # Dams: P1, P2 -> 2 unique dams
  # Parents (union): P1, M1, P2 -> 3 unique parents
  expect_equal(s$n_parents, 3)
})

test_that("selfing: true self-fertilization (Sire == Dam for offspring)", {
  skip_if_not_installed("visPedigree")
  
  # Self-fertilization: C1 has same Sire and Dam (A selfs)
  ped <- data.frame(
    Ind = c("A", "C1", "C2"),
    Sire = c(NA, "A", "A"),
    Dam = c(NA, "A", "A"),
    stringsAsFactors = FALSE
  )
  
  res <- suppressMessages(tidyped(ped, selfing = TRUE))
  
  expect_equal(res[Ind == "A", Sex], "monoecious")
  expect_equal(res[Ind == "C1", Sire], "A")
  expect_equal(res[Ind == "C1", Dam], "A")
  
  # Generation check: A is founder (Gen 1), C1/C2 are offspring (Gen 2)
  expect_equal(res[Ind == "A", Gen], 1L)
  expect_equal(res[Ind == "C1", Gen], 2L)
})

test_that("selfing: inbreeding coefficient computation works", {
  skip_if_not_installed("visPedigree")
  
  # Self-fertilization should produce inbred offspring (f = 0.5)
  ped <- data.frame(
    Ind = c("A", "C1"),
    Sire = c(NA, "A"),
    Dam = c(NA, "A"),
    stringsAsFactors = FALSE
  )
  
  res <- suppressMessages(tidyped(ped, selfing = TRUE, inbreed = TRUE))
  
  expect_equal(res[Ind == "A", f], 0, tolerance = 1e-10)
  expect_equal(res[Ind == "C1", f], 0.5, tolerance = 1e-10)
})

test_that("selfing: multi-generation plant pedigree", {
  skip_if_not_installed("visPedigree")
  
  # Simulate a realistic plant breeding scenario:
  # P1 x P2 -> F1
  # F1 selfs -> F2_1, F2_2
  # F2_1 selfs -> F3
  ped <- data.frame(
    Ind =  c("P1", "P2", "F1", "F2_1", "F2_2", "F3"),
    Sire = c(NA,   NA,   "P1", "F1",   "F1",   "F2_1"),
    Dam =  c(NA,   NA,   "P2", "F1",   "F1",   "F2_1"),
    stringsAsFactors = FALSE
  )
  
  res <- suppressMessages(tidyped(ped, selfing = TRUE, inbreed = TRUE))
  
  # F1 is monoecious (appears as both Sire and Dam)
  expect_equal(res[Ind == "F1", Sex], "monoecious")
  # F2_1 is monoecious (appears as both Sire and Dam for F3)
  expect_equal(res[Ind == "F2_1", Sex], "monoecious")
  # P1 is only sire -> male
  expect_equal(res[Ind == "P1", Sex], "male")
  # P2 is only dam -> female
  expect_equal(res[Ind == "P2", Sex], "female")
  
  # Generation structure
  expect_equal(res[Ind == "P1", Gen], 1L)
  expect_equal(res[Ind == "F1", Gen], 2L)
  expect_equal(res[Ind == "F2_1", Gen], 3L)
  expect_equal(res[Ind == "F3", Gen], 4L)
  
  # Inbreeding: F2_1 and F2_2 from selfing F1, f = 0.5
  expect_equal(res[Ind == "F2_1", f], 0.5, tolerance = 1e-10)
  # F3 from selfing F2_1 (which itself has f=0.5), f = 0.75
  expect_equal(res[Ind == "F3", f], 0.75, tolerance = 1e-10)
})

test_that("selfing: splitped propagates selfing attribute", {
  skip_if_not_installed("visPedigree")
  
  # Two disconnected groups, one with selfing
  ped <- data.frame(
    Ind =  c("A", "C1", "X", "Y", "Z"),
    Sire = c(NA,  "A",  NA,  "X", NA),
    Dam =  c(NA,  "A",  NA,  NA,  NA),
    stringsAsFactors = FALSE
  )
  
  res <- suppressMessages(tidyped(ped, selfing = TRUE))
  expect_true(isTRUE(attr(res, "selfing")))
  
  # splitped should propagate selfing
  groups <- splitped(res)
  
  # Each group should be a valid tidyped
  for (gp in groups) {
    expect_s3_class(gp, "tidyped")
  }
})

test_that("selfing: selfing hint in error message when selfing = FALSE", {
  skip_if_not_installed("visPedigree")
  
  ped <- data.frame(
    Ind = c("A", "C"),
    Sire = c(NA, "A"),
    Dam = c(NA, "A"),
    stringsAsFactors = FALSE
  )
  
  expect_error(
    tidyped(ped, selfing = FALSE),
    "selfing = TRUE"
  )
})
