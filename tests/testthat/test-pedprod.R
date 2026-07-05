test_that("pedprod matches explicit A and Ainv products", {
  tped <- tidyped(small_ped)
  A <- pedmat(tped, method = "A", sparse = FALSE)
  Ainv <- pedmat(tped, method = "Ainv", sparse = FALSE)

  set.seed(20260704)
  x <- stats::rnorm(nrow(tped))
  X <- cbind(First = x, Second = rev(x), Indicator = seq_len(nrow(tped)) %% 2)

  expect_equal(
    unname(pedprod(tped, x, method = "A")),
    unname(drop(A %*% x)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(pedprod(tped, X, method = "A")),
    unname(A %*% X),
    tolerance = 1e-12
  )
  expect_equal(
    unname(pedprod(tped, x, method = "Ainv")),
    unname(drop(Ainv %*% x)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(pedprod(tped, X, method = "Ainv")),
    unname(Ainv %*% X),
    tolerance = 1e-12
  )
})

test_that("pedprod aligns named vectors and matrices and zero-fills omissions", {
  tped <- tidyped(small_ped)
  A <- pedmat(tped, method = "A", sparse = FALSE)

  weights <- c(Z2 = 0.7, A = 0.3)
  full_weights <- setNames(numeric(nrow(tped)), tped$Ind)
  full_weights[names(weights)] <- weights

  result <- pedprod(tped, weights)
  expect_named(result, tped$Ind)
  expect_equal(
    unname(result),
    unname(drop(A %*% full_weights)),
    tolerance = 1e-12
  )

  schemes <- rbind(
    Z2 = c(One = 0.7, Two = 0.2),
    A = c(One = 0.3, Two = 0.8)
  )
  full_schemes <- matrix(
    0,
    nrow = nrow(tped),
    ncol = ncol(schemes),
    dimnames = list(tped$Ind, colnames(schemes))
  )
  full_schemes[rownames(schemes), ] <- schemes

  matrix_result <- pedprod(tped, schemes)
  expect_equal(rownames(matrix_result), tped$Ind)
  expect_equal(colnames(matrix_result), colnames(schemes))
  expect_equal(
    unname(matrix_result),
    unname(A %*% full_schemes),
    tolerance = 1e-12
  )
})

test_that("pedprod handles founders, one known parent, selfing, and disconnected groups", {
  raw <- data.frame(
    Ind = c("F1", "F2", "H1", "S1", "G1", "G2", "H2"),
    Sire = c(NA, NA, "F1", "H1", "S1", NA, "G1"),
    Dam = c(NA, NA, NA, "H1", "F2", NA, "G2"),
    stringsAsFactors = FALSE
  )
  tped <- tidyped(raw, selfing = TRUE)
  A <- pedmat(tped, method = "A", sparse = FALSE)
  Ainv <- pedmat(tped, method = "Ainv", sparse = FALSE)
  x <- seq_len(nrow(tped)) / nrow(tped)

  expect_equal(
    unname(A %*% Ainv),
    diag(nrow(tped)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(pedprod(tped, x)),
    unname(drop(A %*% x)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(pedprod(tped, x, method = "Ainv")),
    unname(drop(Ainv %*% x)),
    tolerance = 1e-12
  )
})

test_that("pedprod accepts raw pedigrees and tidyped objects without numeric columns", {
  raw <- as.data.frame(small_ped)
  tped_without_num <- tidyped(raw, addnum = FALSE)
  x <- rep(1, nrow(tped_without_num))

  raw_result <- pedprod(raw, x)
  tidyped_result <- pedprod(tped_without_num, x)

  expect_equal(raw_result, tidyped_result)
  expect_named(tidyped_result, tidyped(raw)$Ind)
})

test_that("pedprod validates right-hand sides and pedigree structure", {
  tped <- tidyped(small_ped)
  n <- nrow(tped)

  expect_error(pedprod(tped, numeric(n - 1L)), "must have length")
  expect_error(pedprod(tped, matrix(1, nrow = n - 1L)), "must have.*rows")
  expect_error(pedprod(tped, letters[seq_len(n)]), "numeric or logical")
  expect_error(pedprod(tped, replace(numeric(n), 1L, NA_real_)), "finite")
  expect_error(pedprod(tped, replace(numeric(n), 1L, Inf)), "finite")
  expect_error(pedprod(tped, numeric()), "at least one")
  expect_error(pedprod(tped, matrix(numeric(), nrow = n, ncol = 0L)), "at least one")

  unknown <- c(A = 0.5, UNKNOWN = 0.5)
  expect_error(pedprod(tped, unknown), "unknown individual IDs")

  duplicated <- c(A = 0.5, A = 0.5)
  expect_error(pedprod(tped, duplicated), "unique individual IDs")

  missing_name <- c(0.5, 0.5)
  names(missing_name) <- c("A", NA_character_)
  expect_error(pedprod(tped, missing_name), "missing or empty IDs")

  incomplete <- suppressWarnings(tped[-1L])
  expect_error(
    pedprod(incomplete, rep(1, nrow(incomplete))),
    "structurally complete pedigree"
  )

  expect_error(pedprod(splitped(tped), rep(1, n)), "does not support 'splitped'")
  expect_error(pedprod(tped, rep(1, n), method = "D"))
})

test_that("pedprod remains available beyond the dense A size guard", {
  skip_on_cran()

  n <- 25001L
  ids <- paste0("I", seq_len(n))
  raw <- data.frame(
    Ind = ids,
    Sire = c(NA_character_, ids[-n]),
    Dam = NA_character_,
    stringsAsFactors = FALSE
  )
  tped <- suppressMessages(tidyped(raw))

  result <- pedprod(tped, setNames(1, ids[n]))

  expect_length(result, n)
  expect_true(all(is.finite(result)))
  expect_error(pedmat(tped, method = "A"), "too large")
})
