library(testthat)
library(data.table)

test_that("summary.tidyped returns expected core fields", {
  data(simple_ped, package = "visPedigree")
  ped <- tidyped(simple_ped)
  s <- summary(ped)

  expect_s3_class(s, "summary.tidyped")
  expect_equal(s$n_ind, nrow(ped))
  expect_equal(s$n_male + s$n_female + s$n_unknown_sex, s$n_ind)
  expect_true(is.numeric(s$max_gen))
  expect_true(!is.null(s$gen_dist))
  expect_true(is.numeric(s$n_families))
  expect_true(is.numeric(s$n_parents))
})

test_that("summary.tidyped handles addgen=FALSE", {
  data(simple_ped, package = "visPedigree")
  ped <- tidyped(simple_ped, addgen = FALSE)
  s <- summary(ped)

  expect_true(is.na(s$max_gen))
  expect_true(is.na(s$n_isolated))
  expect_null(s$gen_dist)
})

test_that("summary.tidyped handles all-NA inbreeding column", {
  data(simple_ped, package = "visPedigree")
  ped <- tidyped(simple_ped)
  data.table::setalloccol(ped)
  data.table::set(ped, j = "f", value = rep(NA_real_, nrow(ped)))

  s <- summary(ped)
  expect_equal(s$n_inbred, 0)
  expect_true(all(is.na(unlist(s$f_stats))))
})

test_that("summary.tidyped validates input object", {
  bad_class <- data.table(Ind = "A", Sire = NA_character_, Dam = NA_character_, Sex = "male")
  expect_error(summary.tidyped(bad_class), "Object must be of class 'tidyped'")

  bad_cols <- data.table(Ind = "A", Sire = NA_character_, Dam = NA_character_)
  class(bad_cols) <- c("tidyped", class(bad_cols))
  expect_error(summary(bad_cols), "missing required columns: Sex")
})
