test_that("pedhalflife returns correct structure", {
  data(simple_ped)
  tped <- suppressMessages(tidyped(simple_ped))

  res <- suppressMessages(
    suppressWarnings(pedhalflife(tped, timevar = "Gen"))
  )

  # Class and top-level components

  expect_s3_class(res, "pedhalflife")
  expect_true("timeseries" %in% names(res))
  expect_true("decay" %in% names(res))

  # Timeseries columns
  ts <- res$timeseries
  expect_s3_class(ts, "data.table")
  expect_true(all(c("Time", "NRef", "Fe", "Fa", "Fg",
                     "LnFe", "LnFa", "LnFg", "LnFaFe", "LnFgFa",
                     "TimeStep") %in% names(ts)))
  expect_true(nrow(ts) >= 2)

  # Decay columns
  d <- res$decay
  expect_s3_class(d, "data.table")
  expect_true(all(c("LambdaE", "LambdaB", "LambdaD",
                     "LambdaTotal", "THalf") %in% names(d)))
  expect_equal(nrow(d), 1L)
})


test_that("lambda decomposition sums correctly (OLS linearity)", {
  data(simple_ped)
  tped <- suppressMessages(tidyped(simple_ped))

  res <- suppressMessages(
    suppressWarnings(pedhalflife(tped, timevar = "Gen"))
  )
  d <- res$decay

  # By linearity of least squares, lambda_total == lambda_e + lambda_b + lambda_d
  if (!is.na(d$LambdaTotal)) {
    expect_equal(d$LambdaTotal,
                 d$LambdaE + d$LambdaB + d$LambdaD,
                 tolerance = 1e-10)
  }

  # THalf identity
  if (!is.na(d$THalf) && d$LambdaTotal > 0) {
    expect_equal(d$THalf, log(2) / d$LambdaTotal, tolerance = 1e-10)
  }
})


test_that("log columns are consistent with raw values", {
  data(simple_ped)
  tped <- suppressMessages(tidyped(simple_ped))

  res <- suppressMessages(
    suppressWarnings(pedhalflife(tped, timevar = "Gen"))
  )
  ts <- res$timeseries

  expect_equal(ts$LnFe, log(ts$Fe), tolerance = 1e-12)
  expect_equal(ts$LnFa, log(ts$Fa), tolerance = 1e-12)
  expect_equal(ts$LnFg, log(ts$Fg), tolerance = 1e-12)
  expect_equal(ts$LnFaFe, log(ts$Fa / ts$Fe), tolerance = 1e-12)
  expect_equal(ts$LnFgFa, log(ts$Fg / ts$Fa), tolerance = 1e-12)
})


test_that("pedhalflife rejects invalid inputs", {
  data(simple_ped)
  tped <- suppressMessages(tidyped(simple_ped))

  # Missing column
  expect_error(pedhalflife(tped, timevar = "NoSuchCol"),
               "Column 'NoSuchCol' not found")

  # Single time point
  tped2 <- data.table::copy(tped)
  tped2$Year <- 2000
  expect_error(suppressMessages(pedhalflife(tped2, timevar = "Year")),
               "At least two distinct")
})


test_that("print and plot methods work without error", {
  data(simple_ped)
  tped <- suppressMessages(tidyped(simple_ped))

  res <- suppressMessages(
    suppressWarnings(pedhalflife(tped, timevar = "Gen"))
  )

  # print returns invisible x
  out <- capture.output(ret <- print(res))
  expect_identical(ret, res)
  expect_true(length(out) > 0)

  # plot returns a lattice trellis object
  p <- plot(res, type = "log")
  expect_s3_class(p, "trellis")

  p2 <- plot(res, type = "raw")
  expect_s3_class(p2, "trellis")
})
