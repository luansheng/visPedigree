# test-summary-pedmat-vismat-vispstat.R
# Basic tests for summary_pedmat, vismat, and vispstat

# --- summary_pedmat ---
test_that("summary_pedmat works on A matrix", {
  tp <- tidyped(simple_ped)
  A <- pedmat(tp, method = "A")
  s <- summary_pedmat(A)

  expect_s3_class(s, "summary.pedmat")
  expect_equal(s$method, "A")
  expect_false(s$compact)
  expect_equal(s$n_original, nrow(tp))
})

test_that("summary.pedmat S3 dispatch works", {
  tp <- tidyped(simple_ped)
  A <- pedmat(tp, method = "A")
  # Use summary_pedmat directly; generic summary() may dispatch to
  # summary.matrix before summary.pedmat depending on class order.
  s <- summary_pedmat(A)

  expect_s3_class(s, "summary.pedmat")
})

test_that("summary_pedmat errors on non-pedmat input", {
  expect_error(summary_pedmat(matrix(1:4, 2, 2)), "must be a pedmat")
})

# --- vismat ---
test_that("vismat produces a heatmap from pedmat object", {
  tp <- tidyped(simple_ped)
  A <- pedmat(tp, method = "A")

  # vismat returns a lattice/trellis object
  p <- vismat(A, type = "heatmap")
  expect_true(inherits(p, "trellis"))
})

test_that("vismat produces a histogram from pedmat object", {
  tp <- tidyped(simple_ped)
  A <- pedmat(tp, method = "A")

  p <- vismat(A, type = "histogram")
  expect_true(inherits(p, "trellis"))
})

test_that("vismat works with tidyped input", {
  tp <- tidyped(simple_ped)

  p <- vismat(tp, type = "heatmap")
  expect_true(inherits(p, "trellis"))
})

# --- vispstat ---
test_that("vispstat produces ECG histogram", {
  tp <- tidyped(simple_ped)
  stats <- pedstats(tp)

  p <- vispstat(stats, type = "ecg")
  expect_true(inherits(p, "trellis"))
})

test_that("vispstat ECG supports FullGen metric", {
  tp <- tidyped(simple_ped)
  stats <- pedstats(tp)

  p <- vispstat(stats, type = "ecg", metric = "FullGen")
  expect_true(inherits(p, "trellis"))
})

test_that("vispstat errors on non-pedstats input", {
  expect_error(vispstat(list()), "must be a pedstats object")
})

test_that("plot.pedstats dispatches to vispstat", {
  tp <- tidyped(simple_ped)
  stats <- pedstats(tp)

  p <- plot(stats, type = "ecg")
  expect_true(inherits(p, "trellis"))
})
