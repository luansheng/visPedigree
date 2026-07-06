library(testthat)
library(data.table)

# Load package functions if not available (for standalone testing)
if (!exists("tidyped")) {
  # Attempt to locate source files relative to this test file
  # If running from package root: "R/tidyped.R"
  # If running from tests/testthat: "../../R/tidyped.R"
  
  candidates <- c(
    file.path("R", "tidyped.R"),
    file.path("..", "..", "R", "tidyped.R")
  )
  
  for (f in candidates) {
    if (file.exists(f)) {
      # We also need utils for S3 class
      utils_f <- sub("tidyped.R", "utils-s3.R", f)
      inbreed_f <- sub("tidyped.R", "inbreed.R", f) # If inbreed=TRUE tested
      
      if (file.exists(utils_f)) source(utils_f)
      if (file.exists(inbreed_f)) source(inbreed_f)
      source(f)
      break
    }
  }
}

# Load simple_ped data for testing
if (!exists("simple_ped")) {
  data_candidates <- c(
    file.path("data", "simple_ped.rda"),
    file.path("..", "..", "data", "simple_ped.rda")
  )
  for (df in data_candidates) {
    if (file.exists(df)) {
      load(df)
      break
    }
  }
}

test_that("1. Basic Integrity and Types", {
  # Minimal pedigree
  ped_df <- data.frame(
    Ind = c("A", "B", "C"),
    Sire = c(NA, "A", "A"),
    Dam = c(NA, NA, "B")
  )
  
  # Should accept data.frame and return data.table/tidyped
  res <- tidyped(ped_df, addgen = TRUE)
  
  expect_s3_class(res, "tidyped")
  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), 3)
  expect_true("Gen" %in% names(res))
  
  # Check NA handling
  expect_true(is.na(res[Ind == "A", Sire]))
})

test_that("2. Missing Value and Character Normalization", {
  ped_dirty <- data.table(
    Ind = c("X", "Y"),
    Sire = c("0", "X"),    # "0" -> NA
    Dam = c("*", "NA")     # "*" or "NA" -> NA
  )
  
  res <- tidyped(ped_dirty)
  
  # X should have NA parents
  expect_true(is.na(res[Ind == "X", Sire]))
  expect_true(is.na(res[Ind == "X", Dam]))
  
  # Y should have X as Sire
  expect_equal(res[Ind == "Y", Sire], "X")
})

test_that("3. Loop Detection", {
  # A -> B -> A
  ped_loop <- data.frame(
    Ind = c("A", "B"),
    Sire = c("B", NA),
    Dam = c(NA, "A")
  )
  
  expect_error(tidyped(ped_loop), "Pedigree loops detected")
  
  # Self loop A -> A
  ped_self <- data.frame( Ind="A", Sire="A", Dam=NA)
  expect_error(tidyped(ped_self), "Pedigree loops detected")
})

test_that("4. Generation Alignment Logic", {
  # Case: Alignment of introduced founder (E) to mate (C)
  # Gen 1: A, B
  # Gen 2: C (Child of A,B)
  # Founder E mates with C.
  # Gen 3: D (Child of C, E)
  
  # Standard behavior: Founders (no parents) get Gen 1
  # E is a founder, so E gets Gen 1
  # D is child of C (Gen 2) and E (Gen 1), so D gets max(parent_gen) + 1 = Gen 3
  
  ped_align <- data.table(
    Ind = c("A", "B", "C", "D", "E"),
    Sire = c(NA, NA, "A", "C", NA),
    Dam = c(NA, NA, "B", "E", NA)
  )
  
  res <- tidyped(ped_align)
  setkey(res, Ind)
  
  expect_equal(res["C", Gen], 2)
  # E is a founder (no parents), so E gets Gen 1 (standard behavior)
  expect_equal(res["E", Gen], 1)
  # D is child of C (Gen 2) and E (Gen 1) -> Gen 3
  expect_equal(res["D", Gen], 3)
})

test_that("5. Sibling Alignment (Full Sibs)", {
  # Full siblings should have same generation
  # Gen 1: P1, P2
  # Gen 2: C1 (P1xP2)
  # Gen 3: G3 (Child of C1)
  # C2 is full sib of C1 (P1xP2). C2 has no progeny. 
  # Without alignment, C2 might fall to bottom (Gen 3) in bottom-up mode
  # With Sibling/Height logic, C2 matches C1 (Gen 2).
  
  ped_sibs <- data.table(
    Ind = c("P1", "P2", "C1", "C2", "G3"),
    Sire = c(NA, NA, "P1", "P1", "C1"),
    Dam  = c(NA, NA, "P2", "P2", NA)
  )
  
  # Default (top)
  res_top <- tidyped(ped_sibs, genmethod = "top")
  setkey(res_top, Ind)
  expect_equal(res_top["C1", Gen], res_top["C2", Gen], info = "Top mode: Siblings aligned")
  
  # Bottom-up
  res_bottom <- tidyped(ped_sibs, genmethod = "bottom")
  setkey(res_bottom, Ind)
  expect_equal(res_bottom["C1", Gen], 2, info = "Bottom mode: C1 Gen")
  expect_equal(res_bottom["C2", Gen], 2, info = "Bottom mode: C2 aligned with C1")
  expect_equal(res_bottom["G3", Gen], 3, info = "Bottom mode: G3 Gen")
})

test_that("6. Tracing Up/Down", {
  # G1: F1, M1
  # G2: C1 (F1xM1)
  # G3: C2 (C1xUnknown)
  
  ped <- data.table(
    Ind = c("F1", "M1", "C1", "C2"),
    Sire = c(NA, NA, "M1", "C1"),
    Dam = c(NA, NA, "F1", NA)
  )
  
  # Trace UP from C2
  up <- tidyped(ped, cand = "C2", trace = "up")
  expect_true(all(c("F1", "M1", "C1", "C2") %in% up$Ind))
  
  # Trace DOWN from F1
  down <- tidyped(ped, cand = "F1", trace = "down")
  expect_true(all(c("F1", "C1", "C2") %in% down$Ind))
  
  # Trace Gen Limit (e.g., trace up 1 gen from C1)
  # Should include C1 and parents (F1, M1).
  up_lim <- tidyped(ped, cand = "C1", trace = "up", tracegen = 1)
  expect_true(all(c("F1", "M1") %in% up_lim$Ind))
  
  # Trace "all" (Up + Down) from C1
  # Ancestors: F1, M1. Descendants: C2. Self: C1.
  res_all <- tidyped(ped, cand = "C1", trace = "all")
  expect_true(all(c("F1", "M1", "C1", "C2") %in% res_all$Ind))
  expect_equal(sort(res_all$Ind), sort(c("F1", "M1", "C1", "C2")))
  
  # Trace Gen Depth (tracegen > 1) with C2
  # C2 -> C1 -> (F1, M1)
  # tracegen=1 from C2: Should get C2 + Parents (C1). Should missing Grandparents (F1, M1).
  up_gen1 <- tidyped(ped, cand = "C2", trace = "up", tracegen = 1)
  expect_true("C1" %in% up_gen1$Ind)
  expect_false(any(c("F1", "M1") %in% up_gen1$Ind))
  
  # tracegen=2 from C2: Should reach Grandparents
  up_gen2 <- tidyped(ped, cand = "C2", trace = "up", tracegen = 2)
  expect_true(all(c("F1", "M1") %in% up_gen2$Ind))
})

test_that("11. Parameter Validation", {
  ped <- data.table(
    Ind = c("A", "B"),
    Sire = c(NA, NA),
    Dam = c(NA, NA)
  )
  
  # addgen
  expect_error(tidyped(ped, addgen = 1), "The addgen parameter only is assigned using TRUE or FALSE")
  expect_error(tidyped(ped, addgen = "TRUE"), "The addgen parameter only is assigned using TRUE or FALSE")
  
  # addnum
  expect_error(tidyped(ped, addnum = 1), "The addnum parameter only is assigned using TRUE or FALSE")
  
  # inbreed
  expect_error(tidyped(ped, inbreed = "yes"), "The inbreed parameter only is assigned using TRUE or FALSE")
  
  # trace
  expect_error(tidyped(ped, trace = "backwards"), "The trace parameter must be one of: up, down, all")
  
  # tracegen
  expect_error(tidyped(ped, tracegen = "2"), "The tracegen parameter must be a single numeric value")
  
  # ped validation
  expect_error(tidyped(), "The ped parameter cannot be NULL or missing")
  expect_error(tidyped(NULL), "The ped parameter cannot be NULL or missing")
  expect_error(tidyped(data.frame()), "The ped parameter cannot be empty")
})

test_that("10. Deep Pedigree Tests (simple_ped)", {
  skip_if_not(exists("simple_ped"), "simple_ped data not available")
  
  # Ensure simple_ped is a data.table or data.frame
  ped_deep <- as.data.table(simple_ped)
  
  # Run full processing
  res <- tidyped(ped_deep, addgen = TRUE, addnum = TRUE)
  
  # 1. Structure check
  expect_s3_class(res, "tidyped")
  expect_true("Gen" %in% names(res))
  
  # 2. Generation Depth Check
  # Based on previous context, simple_ped has around 6 generations
  max_gen <- max(res$Gen)
  expect_gt(max_gen, 4) 
  
  # 3. Parent-Offspring Generation Constraint
  # Gen(Child) should generally be > Gen(Parent)
  # (Though with overlap generations and heurisitcs, strictly > is expected)
  
  # Merge Sire Gen
  res_sire <- merge(res[, .(IndNum, Gen, SireNum)], res[, .(IndNum, SireGen=Gen)], 
                    by.x="SireNum", by.y="IndNum", all.x=TRUE) # Use Num for speed/safety
  # Compare only known parents
  valid_sire <- res_sire[SireNum > 0]
  
  # Check if Child Gen > Sire Gen
  # Note: Due to alignment heuristics, Child Gen could be much larger, but never smaller or equal?
  # Actually, graph is DAG, so Child Gen MUST be > Sire Gen if purely topological.
  # The heuristic might pull parents UP (increase numeric value) -> making parent older.
  # Small Gen number = Old (Founders). Large Gen number = Young.
  # WAIT: Previous heuristic implementation:
  # Founders aligned to mate's generation (min mate gen).
  # Height propagation: Gen = MaxH - H + 1.
  # Topo sort implies Parents appear before Children.
  # If Gen 1 is old, Gen 6 is young.
  # Then Gen(Child) > Gen(Parent).
  
  expect_true(all(valid_sire$Gen > valid_sire$SireGen))
  
  # 4. Tracing a deep individual
  # Using J5X804 (Generation 6 in prior context)
  cand_id <- "J5X804"
  if (cand_id %in% res$Ind) {
    trace_res <- tidyped(ped_deep, cand = cand_id, trace = "up")
    
    # Should contain cand
    expect_true(cand_id %in% trace_res$Ind)
    # Should contain ancestors (at least one parent if known)
    # J5X804 has parents J4Y326/J4E185 in context
    parents <- c("J4Y326", "J4E185")
    if (any(parents %in% ped_deep$Ind)) {
       expect_true(any(parents %in% trace_res$Ind))
    }
  }
})

test_that("12. addgen=FALSE functionality", {
  # Regression test for bug where addgen=FALSE returned empty result
  ped <- data.table(
    Ind = c("C", "A", "B"),
    Sire = c("A", NA, NA),
    Dam = c("B", NA, NA)
  )
  
  res <- tidyped(ped, addgen = FALSE)
  
  expect_equal(nrow(res), 3)
  expect_false("Gen" %in% names(res))
  
  # Check topological sorting
  # Parents (A, B) should appear before Child (C)
  # But exact order of A vs B is not guaranteed, just before C.
  idx_A <- which(res$Ind == "A")
  idx_B <- which(res$Ind == "B")
  idx_C <- which(res$Ind == "C")
  
  expect_lt(idx_A, idx_C)
  expect_lt(idx_B, idx_C)
})
