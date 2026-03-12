# Test edge cases and boundary conditions

test_that("pedmat rejects multiple methods", {
  skip_if_not_installed("visPedigree")
  tped <- tidyped(small_ped)
  
  expect_error(
    pedmat(tped, method = c("A", "D")),
    "Only a single method may be requested"
  )
})

test_that("pedmat rejects invalid methods", {
  skip_if_not_installed("visPedigree")
  tped <- tidyped(small_ped)
  
  expect_error(
    pedmat(tped, method = "invalid"),
    "Invalid method"
  )
})

test_that("query_relationship rejects inverse matrices", {
  skip_if_not_installed("visPedigree")
  tped <- tidyped(small_ped)
  Ainv <- pedmat(tped, method = "Ainv")
  
  expect_error(
    query_relationship(Ainv, "A", "B"),
    "does not support inverse matrices"
  )
})

test_that("vismat rejects inverse matrices", {
  skip_if_not_installed("visPedigree")
  tped <- tidyped(small_ped)
  Ainv <- pedmat(tped, method = "Ainv")
  
  expect_error(
    vismat(Ainv),
    "does not support inverse matrices"
  )
})

test_that("compact mode expand_pedmat works for vectors", {
  skip_if_not_installed("visPedigree")
  tped <- tidyped(small_ped)
  
  # Get inbreeding coefficients in compact mode
  f_compact <- pedmat(tped, method = "f", compact = TRUE)
  f_full <- pedmat(tped, method = "f", compact = FALSE)
  
  # Expand compact version
  f_expanded <- expand_pedmat(f_compact)
  
  # Should have same length as original pedigree
  expect_equal(length(f_expanded), nrow(tped))
  
  # Names should match original individuals
  expect_equal(names(f_expanded), tped$Ind)
})

test_that("compact mode query_relationship handles full-siblings", {
  skip_if_not_installed("visPedigree")
  tped <- tidyped(small_ped)
  
  # C, D, E are full-siblings (AxB)
  A_compact <- pedmat(tped, method = "A", compact = TRUE)
  A_full <- pedmat(tped, method = "A", compact = FALSE)
  
  # Full-sibling relationship should be 0.5
  expect_equal(query_relationship(A_compact, "C", "D"), 0.5)
  expect_equal(query_relationship(A_compact, "C", "E"), 0.5)
  expect_equal(query_relationship(A_compact, "D", "E"), 0.5)
  
  # Should match full matrix
  expect_equal(
    query_relationship(A_compact, "C", "D"),
    query_relationship(A_full, "C", "D")
  )
})

test_that("D matrix compact mode full-sibling query works", {
  skip_if_not_installed("visPedigree")
  tped <- tidyped(small_ped)
  
  # D matrix in compact mode
  D_compact <- pedmat(tped, method = "D", compact = TRUE)
  D_full <- pedmat(tped, method = "D", compact = FALSE)
  
  # Query full-siblings
  d_cd_compact <- query_relationship(D_compact, "C", "D")
  d_cd_full <- query_relationship(D_full, "C", "D")
  
  # Should match (or be close due to formula)
  expect_equal(d_cd_compact, d_cd_full, tolerance = 1e-6)
})

test_that("sex conflict detection works", {
  skip_if_not_installed("visPedigree")
  
  # Create pedigree with sex conflict
  bad_ped <- data.frame(
    Ind = c("A", "B", "C"),
    Sire = c(NA, NA, "A"),
    Dam = c(NA, NA, "A"),  # A appears as both Sire and Dam!
    stringsAsFactors = FALSE
  )
  
  expect_error(
    tidyped(bad_ped),
    "Sex conflict detected"
  )
  
  # Same pedigree should work with selfing = TRUE
  res <- tidyped(bad_ped, selfing = TRUE)
  expect_s3_class(res, "tidyped")
  expect_equal(res[Ind == "A", Sex], "monoecious")
  expect_true(isTRUE(attr(res, "selfing")))
})

test_that("sex annotation conflict detection works", {
  skip_if_not_installed("visPedigree")
  
  # Create pedigree with explicit sex conflict
  bad_ped <- data.frame(
    Ind = c("A", "B", "C"),
    Sire = c(NA, NA, "A"),
    Dam = c(NA, NA, "B"),
    Sex = c("female", "male", NA),  # A is female but used as Sire!
    stringsAsFactors = FALSE
  )
  
  expect_error(
    tidyped(bad_ped),
    "Sex annotation conflicts"
  )
})

test_that("query_relationship id2=NULL returns correctly named vector", {
  skip_if_not_installed("visPedigree")
  tped <- tidyped(small_ped)
  A_compact <- pedmat(tped, method = "A", compact = TRUE)
  
  # Query row for individual A
  row_a <- query_relationship(A_compact, "A")
  
  # Should have names
  expect_true(!is.null(names(row_a)))
  
  # Length should match matrix dimension
  mat <- A_compact
  class(mat) <- setdiff(class(mat), "pedmat")
  expect_equal(length(row_a), ncol(mat))
  
  # Names should match matrix column names
  expect_equal(names(row_a), colnames(mat))
})
