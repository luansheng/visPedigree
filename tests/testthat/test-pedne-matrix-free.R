test_that("matrix-free sampled statistics match an explicit relationship matrix", {
  ped <- suppressMessages(tidyped(simple_ped))
  ecg <- pedecg(ped)
  ped <- merge(ped, ecg[, .(Ind, ECG)], by = "Ind", all.x = TRUE)
  data.table::setorder(ped, IndNum)

  target_idx <- which(ped$Gen == 2L)
  f_res <- cpp_calculate_inbreeding(ped$SireNum, ped$DamNum)
  actual <- calc_sampled_relationship_stats(
    sire = ped$SireNum,
    dam = ped$DamNum,
    dii = f_res$dii,
    f = f_res$f,
    target_idx = target_idx,
    ecg = ped$ECG,
    need_delta = TRUE,
    need_mean = TRUE,
    batch_size = 3L
  )

  A <- as.matrix(pedmat(ped, method = "A", sparse = FALSE))
  A_target <- A[target_idx, target_idx, drop = FALSE]
  pair_mask <- lower.tri(A_target)
  pair_ecg <- outer(ped$ECG[target_idx], ped$ECG[target_idx], "+") / 2
  valid <- pair_mask & pair_ecg > 0
  expected_delta <- mean(
    1 - (1 - A_target[valid] / 2)^(1 / pair_ecg[valid])
  )

  expect_equal(actual$DeltaC, expected_delta, tolerance = 1e-12)
  expect_equal(actual$MeanRelOff, mean(A_target[pair_mask]), tolerance = 1e-12)
  expect_equal(actual$MeanF, mean(diag(A_target) - 1), tolerance = 1e-12)
})

test_that("coancestry metric selection avoids unnecessary calculations", {
  ped <- suppressMessages(tidyped(simple_ped))
  ecg <- pedecg(ped)
  ped <- merge(ped, ecg[, .(Ind, ECG)], by = "Ind", all.x = TRUE)
  data.table::setorder(ped, IndNum)
  ped$Cohort <- ped$Gen

  ne_only <- suppressWarnings(calc_ne_coancestry(
    ped, ped, "Cohort", nsamples = 1000, seed = 42L, metrics = "ne"
  ))
  fg_only <- suppressWarnings(calc_ne_coancestry(
    ped, ped, "Cohort", nsamples = 1000, seed = 42L, metrics = "fg"
  ))
  both <- suppressWarnings(calc_ne_coancestry(
    ped, ped, "Cohort", nsamples = 1000, seed = 42L,
    metrics = c("ne", "fg")
  ))

  expect_equal(ne_only$DeltaC, both$DeltaC, tolerance = 1e-12)
  expect_equal(ne_only$Ne, both$Ne, tolerance = 1e-12)
  expect_true(all(is.na(ne_only$MeanCoan)))
  expect_true(all(is.na(ne_only$fg)))
  expect_true(all(is.na(ne_only$NSampledCoan)))

  expect_equal(fg_only$MeanCoan, both$MeanCoan, tolerance = 1e-12)
  expect_equal(fg_only$fg, both$fg, tolerance = 1e-12)
  expect_true(all(is.na(fg_only$DeltaC)))
  expect_true(all(is.na(fg_only$Ne)))
})

test_that("single sampled individual preserves missing pair statistics", {
  ped <- suppressMessages(tidyped(simple_ped))
  f_res <- cpp_calculate_inbreeding(ped$SireNum, ped$DamNum)
  target_idx <- which.max(ped$Gen)

  actual <- calc_sampled_relationship_stats(
    sire = ped$SireNum,
    dam = ped$DamNum,
    dii = f_res$dii,
    f = f_res$f,
    target_idx = target_idx,
    ecg = rep(1, nrow(ped))
  )

  expect_true(is.na(actual$DeltaC))
  expect_true(is.na(actual$MeanRelOff))
  expect_equal(actual$MeanF, f_res$f[target_idx])
})
