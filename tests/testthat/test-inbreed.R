test_that("inbreed works correctly", {
  ped <- data.frame(
    Ind = c("A", "B", "C", "D", "E", "F"),
    Sire = c(NA, NA, "A", "A", "C", "C"),
    Dam = c(NA, NA, "B", "B", "D", "D")
  )
  tped <- tidyped(ped)
  
  res <- inbreed(tped)
  
  expect_s3_class(res, "tidyped")
  expect_true("f" %in% names(res))
  
  # A and B are founders, f = 0
  expect_equal(res$f[res$Ind == "A"], 0)
  expect_equal(res$f[res$Ind == "B"], 0)
  
  # C and D are offspring of A and B, f = 0
  expect_equal(res$f[res$Ind == "C"], 0)
  expect_equal(res$f[res$Ind == "D"], 0)
  
  # E and F are offspring of full siblings C and D, f = 0.25
  expect_equal(res$f[res$Ind == "E"], 0.25)
  expect_equal(res$f[res$Ind == "F"], 0.25)
})

test_that("inbreed handles missing parents", {
  ped <- data.frame(
    Ind = c("A", "B", "C", "D"),
    Sire = c(NA, NA, "A", "C"),
    Dam = c(NA, NA, NA, "B")
  )
  tped <- tidyped(ped)
  
  res <- inbreed(tped)
  
  expect_s3_class(res, "tidyped")
  expect_true("f" %in% names(res))
  expect_equal(res$f, c(0, 0, 0, 0))
})

test_that("inbreed validates input", {
  # ensure_tidyped errors on data.frame missing core columns
  expect_error(inbreed(data.frame(Ind="A", Sire=NA, Dam=NA)), "tidyped object")
})
