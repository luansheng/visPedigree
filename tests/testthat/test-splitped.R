test_that("splitped works correctly on connected pedigree", {
  ped <- data.frame(
    Ind = c("A", "B", "C", "D"),
    Sire = c(NA, NA, "A", "A"),
    Dam = c(NA, NA, "B", "B")
  )
  tped <- tidyped(ped)
  
  res <- splitped(tped)
  
  expect_s3_class(res, "splitped")
  expect_equal(attr(res, "n_groups"), 1)
  expect_equal(attr(res, "total"), 4)
  expect_equal(attr(res, "n_isolated"), 0)
  expect_equal(length(res), 1)
  expect_s3_class(res$GP1, "tidyped")
  expect_equal(nrow(res$GP1), 4)
})

test_that("splitped works correctly on disconnected pedigree", {
  ped <- data.frame(
    Ind = c("A", "B", "C", "D", "E", "F", "G"),
    Sire = c(NA, NA, "A", "A", NA, NA, "E"),
    Dam = c(NA, NA, "B", "B", NA, NA, "F")
  )
  tped <- tidyped(ped)
  
  res <- splitped(tped)
  
  expect_s3_class(res, "splitped")
  expect_equal(attr(res, "n_groups"), 2)
  expect_equal(attr(res, "total"), 7)
  expect_equal(attr(res, "n_isolated"), 0)
  expect_equal(length(res), 2)
  
  # Groups are sorted by size descending
  expect_equal(nrow(res$GP1), 4)
  expect_equal(nrow(res$GP2), 3)
  
  expect_true(all(res$GP1$Ind %in% c("A", "B", "C", "D")))
  expect_true(all(res$GP2$Ind %in% c("E", "F", "G")))
})

test_that("splitped handles isolated individuals", {
  ped <- data.frame(
    Ind = c("A", "B", "C", "D", "E", "F"),
    Sire = c(NA, NA, "A", NA, NA, NA),
    Dam = c(NA, NA, "B", NA, NA, NA)
  )
  tped <- tidyped(ped)
  
  res <- splitped(tped)
  
  expect_s3_class(res, "splitped")
  expect_equal(attr(res, "n_groups"), 1)
  expect_equal(attr(res, "total"), 3)
  expect_equal(attr(res, "n_isolated"), 3)
  expect_equal(sort(attr(res, "isolated")), c("D", "E", "F"))
  expect_equal(nrow(res$GP1), 3)
})

test_that("splitped handles all isolated individuals", {
  # Cannot construct an empty tidyped with tidyped() because "All parents are missing"
  # So we manually construct a mimic tidyped object for testing splitped() logic
  
  ped <- data.table::data.table(
    Ind = c("A", "B", "C"),
    Sire = rep(NA_character_, 3),
    Dam = rep(NA_character_, 3),
    Sex = rep(NA_character_, 3),
    Gen = rep(0L, 3), # Isolated individuals have Gen=0
    IndNum = 1:3,
    SireNum = rep(0L, 3),
    DamNum = rep(0L, 3)
  )
  class(ped) <- c("tidyped", "data.table", "data.frame")
  
  res <- splitped(ped)
  
  expect_s3_class(res, "splitped")
  expect_equal(attr(res, "n_groups"), 0)
  expect_equal(attr(res, "total"), 0)
  expect_equal(attr(res, "n_isolated"), 3)
  expect_equal(length(res), 0)
})

test_that("splitped input validation", {
  expect_error(splitped(data.frame(Ind="A", Sire=NA, Dam=NA)), "tidyped object")
})

test_that("print and summary methods for splitped", {
  ped <- data.frame(
    Ind = c("A", "B", "C", "D", "E", "F", "G", "H"),
    Sire = c(NA, NA, "A", "A", NA, NA, "E", NA),
    Dam = c(NA, NA, "B", "B", NA, NA, "F", NA)
  )
  tped <- tidyped(ped)
  res <- splitped(tped)
  
  expect_output(print(res), "Pedigree Split Result")
  expect_output(print(res), "Total individuals in groups: 7")
  expect_output(print(res), "Isolated individuals \\(Gen=0\\): 1")
  expect_output(print(res), "Number of groups:  2")
  
  expect_output(summary(res), "Summary of Pedigree Split")
  expect_output(summary(res), "Grand total: 8")
  expect_output(summary(res), "Connectivity: Pedigree contains disconnected groups")
})
