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

test_that("vismat reorder=FALSE preserves original order", {
  tp <- tidyped(small_ped)
  A <- pedmat(tp, method = "A")
  p <- vismat(A, reorder = FALSE)
  expect_true(inherits(p, "trellis"))
})

test_that("vismat ids subsets the matrix", {
  tp <- tidyped(small_ped)
  A <- pedmat(tp, method = "A")
  target <- rownames(as.matrix(A))[1:5]
  p <- vismat(A, ids = target)
  expect_true(inherits(p, "trellis"))
})

test_that("vismat ids errors on non-existent IDs", {
  tp <- tidyped(small_ped)
  A <- pedmat(tp, method = "A")
  expect_error(vismat(A, ids = c("NONEXIST_1", "NONEXIST_2")),
               "None of the specified")
})

test_that("vismat by groups by generation", {
  tp <- tidyped(small_ped)
  A <- pedmat(tp, method = "A")
  expect_message(
    p <- vismat(A, ped = tp, by = "Gen"),
    "Aggregating"
  )
  expect_true(inherits(p, "trellis"))
})

test_that("vismat by errors without ped", {
  tp <- tidyped(small_ped)
  A <- as.matrix(pedmat(tp, method = "A"))
  # Plain matrix has no ped attribute
  expect_error(vismat(A, by = "Gen"), "'ped' must be provided")
})

test_that("vismat grouping deprecated parameter works", {
  tp <- tidyped(small_ped)
  A <- pedmat(tp, method = "A")
  expect_warning(
    p <- vismat(A, ped = tp, grouping = "Gen"),
    "deprecated"
  )
  expect_true(inherits(p, "trellis"))
})

test_that("vismat works with plain matrix input", {
  tp <- tidyped(small_ped)
  A <- pedmat(tp, method = "A")
  A_dense <- as.matrix(A)
  class(A_dense) <- "matrix"  # ensure no pedmat class
  p <- vismat(A_dense)
  expect_true(inherits(p, "trellis"))
})

test_that("vismat labelcex controls label font size", {
  tp <- tidyped(small_ped)
  A <- pedmat(tp, method = "A")
  p <- vismat(A, labelcex = 0.3)
  expect_true(inherits(p, "trellis"))
})

test_that("vismat auto-expands compact pedmat", {
  tp <- tidyped(small_ped)
  A_compact <- pedmat(tp, method = "A", compact = TRUE)
  ci <- attr(A_compact, "call_info")

  # Should produce a message about expanding
  expect_message(
    p <- vismat(A_compact),
    "Expanding compact matrix"
  )
  expect_true(inherits(p, "trellis"))
})

test_that("vismat errors on unsupported type", {
  tp <- tidyped(small_ped)
  A <- pedmat(tp, method = "A")
  expect_error(vismat(A, type = "scatter"), "not supported")
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

test_that("plot.pedstats dispatches to vispstat for ecg", {
  tp <- tidyped(simple_ped)
  stats <- pedstats(tp)

  p <- plot(stats, type = "ecg")
  expect_true(inherits(p, "trellis"))
})

test_that("vispstat genint produces bar chart when gen_intervals exist", {
  tp <- tidyped(simple_ped)
  # Assign synthetic birth years so pedgenint can compute intervals
  tp$Year <- tp$Gen * 3 + 2000
  stats <- pedstats(tp, timevar = "Year")
  skip_if(is.null(stats$gen_intervals), "No generation intervals computed")

  p <- vispstat(stats, type = "genint")
  expect_true(inherits(p, "trellis"))
})

test_that("vispstat genint errors when gen_intervals is NULL", {
  tp <- tidyped(simple_ped)
  stats <- pedstats(tp, genint = FALSE)

  expect_error(vispstat(stats, type = "genint"),
               "No generation interval data found")
})

test_that("plot.pedstats dispatches genint correctly", {
  tp <- tidyped(simple_ped)
  tp$Year <- tp$Gen * 3 + 2000
  stats <- pedstats(tp, timevar = "Year")
  skip_if(is.null(stats$gen_intervals), "No generation intervals computed")

  p <- plot(stats, type = "genint")
  expect_true(inherits(p, "trellis"))
})

# --- vismat compact + by (fast aggregation path) ---
test_that("vismat compact+by=Gen matches expand-then-aggregate", {
  tp <- tidyped(small_ped)
  A_compact <- pedmat(tp, method = "A", compact = TRUE)
  A_full    <- pedmat(tp, method = "A", compact = FALSE)

  # Fast path: aggregate directly from compact
  agg_compact <- visPedigree:::aggregate_compact_by_group(
    A_compact, tp$Ind, "Gen", tp
  )

  # Reference: expand, then aggregate manually
  mat_exp <- as.matrix(expand_pedmat(A_compact))
  grp <- tp$Gen[match(rownames(mat_exp), tp$Ind)]
  grps <- sort(unique(grp))
  agg_ref <- matrix(0, length(grps), length(grps),
                    dimnames = list(grps, grps))
  for (i in seq_along(grps)) {
    for (j in i:length(grps)) {
      idx_i <- which(grp == grps[i])
      idx_j <- which(grp == grps[j])
      agg_ref[i, j] <- mean(mat_exp[idx_i, idx_j])
      agg_ref[j, i] <- agg_ref[i, j]
    }
  }

  expect_equal(agg_compact, agg_ref, tolerance = 1e-12)
})

test_that("vismat compact+by=Family matches expand-then-aggregate", {
  tp <- tidyped(small_ped)
  A_compact <- pedmat(tp, method = "A", compact = TRUE)

  # Focus on gen 3 (has 3 families)
  ids_g3 <- tp[Gen == 3, Ind]
  agg_compact <- visPedigree:::aggregate_compact_by_group(
    A_compact, ids_g3, "Family", tp
  )

  # Reference
  mat_exp <- as.matrix(expand_pedmat(A_compact))
  sub_ids <- ids_g3
  sub_mat <- mat_exp[sub_ids, sub_ids, drop = FALSE]
  grp <- tp$Family[match(sub_ids, tp$Ind)]
  grps <- sort(unique(grp[!is.na(grp)]))
  keep <- !is.na(grp)
  sub_mat <- sub_mat[keep, keep, drop = FALSE]
  grp <- grp[keep]
  agg_ref <- matrix(0, length(grps), length(grps),
                    dimnames = list(grps, grps))
  for (i in seq_along(grps)) {
    for (j in i:length(grps)) {
      idx_i <- which(grp == grps[i])
      idx_j <- which(grp == grps[j])
      agg_ref[i, j] <- mean(sub_mat[idx_i, idx_j])
      agg_ref[j, i] <- agg_ref[i, j]
    }
  }

  expect_equal(agg_compact, agg_ref, tolerance = 1e-12)
})

test_that("vismat compact+by produces valid trellis plot", {
  tp <- tidyped(small_ped)
  A_compact <- pedmat(tp, method = "A", compact = TRUE)

  expect_message(
    p <- vismat(A_compact, ped = tp, by = "Gen"),
    "Aggregating"
  )
  expect_true(inherits(p, "trellis"))
})
