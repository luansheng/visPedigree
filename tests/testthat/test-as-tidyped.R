# test-as-tidyped.R
# Tests for as_tidyped(), ensure_tidyped(), and S3 class recovery

test_that("as_tidyped restores class after rbind", {
  tp <- tidyped(simple_ped)
  tp2 <- rbind(tp[1:5], tp[6:10])

  expect_false(inherits(tp2, "tidyped"))
  expect_true(inherits(tp2, "data.table"))

  tp3 <- as_tidyped(tp2)
  expect_s3_class(tp3, "tidyped")
  expect_s3_class(tp3, "data.table")
})

test_that("as_tidyped restores class from plain data.frame", {
  tp <- tidyped(simple_ped)
  tp_df <- as.data.frame(tp)
  tp_sub <- tp_df[tp_df$Gen > 1, ]

  expect_false(inherits(tp_sub, "tidyped"))
  expect_false(inherits(tp_sub, "data.table"))

  tp_restored <- as_tidyped(tp_sub)
  expect_s3_class(tp_restored, "tidyped")
  expect_s3_class(tp_restored, "data.table")
})

test_that("as_tidyped rebuilds IndNum after rbind with overlapping indices", {

  tp <- tidyped(simple_ped)
  # rbind two slices that have overlapping IndNum values
  tp2 <- rbind(tp[1:5], tp[6:10])
  tp3 <- as_tidyped(tp2)

  # IndNum must be sequential 1:nrow

  expect_equal(tp3$IndNum, seq_len(nrow(tp3)))

  # SireNum/DamNum must reference valid IndNum or be 0
  valid_range <- c(0L, tp3$IndNum)
  expect_true(all(tp3$SireNum %in% valid_range))
  expect_true(all(tp3$DamNum %in% valid_range))
})

test_that("as_tidyped rebuilds IndNum from data.frame subset", {
  tp <- tidyped(simple_ped)
  tp_df <- as.data.frame(tp)
  tp_sub <- tp_df[tp_df$Gen > 1, ]
  tp_restored <- as_tidyped(tp_sub)

  expect_equal(tp_restored$IndNum, seq_len(nrow(tp_restored)))
  expect_true(all(tp_restored$SireNum >= 0L))
  expect_true(all(tp_restored$DamNum >= 0L))

  # SireNum should correctly map to Ind in restored object
  for (i in seq_len(nrow(tp_restored))) {
    if (tp_restored$SireNum[i] > 0L) {
      expect_equal(
        tp_restored$Ind[tp_restored$SireNum[i]],
        tp_restored$Sire[i]
      )
    }
    if (tp_restored$DamNum[i] > 0L) {
      expect_equal(
        tp_restored$Ind[tp_restored$DamNum[i]],
        tp_restored$Dam[i]
      )
    }
  }
})

test_that("as_tidyped is a no-op on valid tidyped", {
  tp <- tidyped(simple_ped)
  tp2 <- as_tidyped(tp)

  expect_identical(tp, tp2)
})

test_that("as_tidyped errors on missing core columns", {
  df <- data.frame(Ind = c("A", "B"), Sire = c(NA, "A"))

  expect_error(as_tidyped(df), "Missing columns")
})

test_that("as_tidyped errors on non-data.frame input", {
  expect_error(as_tidyped(1:10), "not a data.frame")
  expect_error(as_tidyped("abc"), "not a data.frame")
})

test_that("ensure_tidyped auto-recovers with message", {
  tp <- tidyped(simple_ped)
  tp2 <- rbind(tp[1:5], tp[6:10])

  expect_message(
    tp3 <- ensure_tidyped(tp2),
    "Restoring automatically"
  )
  expect_s3_class(tp3, "tidyped")
  expect_equal(tp3$IndNum, seq_len(nrow(tp3)))
})

test_that("ensure_tidyped is silent on valid tidyped", {
  tp <- tidyped(simple_ped)
  expect_silent(tp2 <- ensure_tidyped(tp))
  expect_identical(tp, tp2)
})

test_that("ensure_tidyped errors on non-data.frame", {
  expect_error(ensure_tidyped(1:10), "tidyped object")
})

test_that("ensure_tidyped errors on missing core columns", {
  df <- data.frame(Ind = c("A", "B"), Sire = c(NA, "A"))
  expect_error(ensure_tidyped(df), "tidyped object")
})

test_that("analysis functions work after class loss via rbind", {
  tp <- tidyped(simple_ped)
  tp2 <- rbind(tp[1:10], tp[11:20])

  # pedstats should auto-recover and work
  expect_message(
    stats <- pedstats(tp2),
    "Restoring automatically"
  )
  expect_s3_class(stats, "pedstats")
  expect_true(is_tidyped(tp2))

  # pedecg should then work on the already-restored object
  expect_silent(ecg <- pedecg(tp2))
  expect_true(is.data.table(ecg))
  expect_true("ECG" %in% names(ecg))
})

test_that("splitped works after class loss", {
  tp <- tidyped(small_ped)
  tp2 <- data.table::as.data.table(tp)
  class(tp2) <- c("data.table", "data.frame")  # strip tidyped class

  expect_message(
    result <- splitped(tp2),
    "Restoring automatically"
  )
  expect_s3_class(result, "splitped")
})

test_that("inbreed works after class loss", {
  tp <- tidyped(simple_ped)
  tp2 <- data.table::as.data.table(tp)
  class(tp2) <- c("data.table", "data.frame")  # strip tidyped class

  expect_message(
    result <- inbreed(tp2),
    "Restoring automatically"
  )
  expect_s3_class(result, "tidyped")
  expect_true("f" %in% names(result))
})
