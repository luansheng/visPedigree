library(testthat)
library(data.table)

test_that("[.tidyped degrades incomplete row subsets with warning", {
  ped <- data.table(
    Ind = c("A", "B", "C", "D"),
    Sire = c(NA, NA, "A", "C"),
    Dam = c(NA, NA, "B", "B"),
    Year = c(2000, 2000, 2005, 2006)
  )

  tp <- tidyped(ped)

  expect_warning(
    sub_dt <- tp[Year > 2005],
    "Subsetting removed parent records"
  )

  expect_identical(class(sub_dt), c("data.table", "data.frame"))
  expect_false(is_tidyped(sub_dt))
})

test_that("[.tidyped preserves complete pedigree subsets and rebuilds indices", {
  ped <- data.table(
    Ind = c("A", "B", "C", "D", "E", "F"),
    Sire = c(NA, NA, "A", NA, NA, "D"),
    Dam = c(NA, NA, "B", NA, NA, "E")
  )

  tp <- tidyped(ped)
  sub_tp <- tp[Ind %in% c("A", "B", "C")]

  expect_true(is_tidyped(sub_tp))
  expect_identical(sub_tp$IndNum, 1:3)
  expect_equal(sub_tp[Ind == "C", SireNum], sub_tp[Ind == "A", IndNum])
  expect_equal(sub_tp[Ind == "C", DamNum], sub_tp[Ind == "B", IndNum])
})

test_that(":= keeps tidyped class and metadata", {
  tp <- tidyped(simple_ped)

  tp[, phenotype := seq_len(.N)]

  expect_true(is_tidyped(tp))
  expect_true("phenotype" %in% names(tp))
  expect_identical(pedmeta(tp)$genmethod, "top")
})

test_that("tidyped fast path matches full tracing result", {
  tp_master <- tidyped(simple_ped)

  res_raw <- tidyped(simple_ped, cand = "J5X804", trace = "up", tracegen = 2)
  res_fast <- tidyped(tp_master, cand = "J5X804", trace = "up", tracegen = 2)

  expect_true(is_tidyped(res_fast))
  expect_true(has_candidates(res_fast))
  expect_equal(as.list(res_fast), as.list(res_raw))
  expect_equal(pedmeta(res_fast), pedmeta(res_raw))
})

test_that("state accessors report tidyped contents correctly", {
  tp <- tidyped(simple_ped)
  tp_f <- inbreed(tp)
  tp_c <- tidyped(tp, cand = "J5X804", trace = "up")

  expect_true(is_tidyped(tp))
  expect_false(has_inbreeding(tp))
  expect_true(has_inbreeding(tp_f))
  expect_false(has_candidates(tp))
  expect_true(has_candidates(tp_c))
})

test_that("splitped returns objects while pedsubpop returns summary table", {
  ped <- data.table(
    Ind = c("A", "B", "C", "D", "E", "F", "G"),
    Sire = c(NA, NA, "A", "A", NA, NA, "E"),
    Dam = c(NA, NA, "B", "B", NA, NA, "F")
  )

  tp <- tidyped(ped)
  splits <- splitped(tp)
  stats <- pedsubpop(tp)

  expect_s3_class(splits, "splitped")
  expect_true(all(vapply(splits, is_tidyped, logical(1))))
  expect_s3_class(stats, "data.table")
  expect_false("subpop" %in% names(stats))
  expect_true(all(c("Group", "N", "N_Sire", "N_Dam", "N_Founder") %in% names(stats)))
})

test_that("completeness-sensitive analyses error on truncated subsets", {
  tp <- tidyped(simple_ped)

  expect_warning(
    tp_sub <- tp[Gen > 2],
    "Subsetting removed parent records"
  )

  expect_false(is_tidyped(tp_sub))

  expect_error(
    inbreed(tp_sub),
    "structurally complete pedigree"
  )

  expect_error(
    pedecg(tp_sub),
    "structurally complete pedigree"
  )

  expect_error(
    pedmat(tp_sub, method = "f"),
    "structurally complete pedigree"
  )
})