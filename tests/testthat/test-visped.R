library(testthat)
library(data.table)
library(visPedigree)

test_that("visped works with basic tidy input", {
  tidy_ped <- tidyped(simple_ped, addgen = TRUE, addnum = TRUE)
  # Check return structure
  res <- visped(tidy_ped, showgraph = FALSE, file = tempfile())
  expect_type(res, "list")
  expect_named(res, c("g", "layout", "canvas_width", "canvas_height", "node_size", "best_cex"))
  expect_s3_class(res$g, "igraph")
})

test_that("visped handles raw input by auto-tidying", {
  # Raw data frame
  res <- visped(simple_ped, showgraph = FALSE, file = tempfile())
  expect_s3_class(res$g, "igraph")
  
  # Check dimensions - simple_ped has 52 individuals in example
  expect_gt(igraph::vcount(res$g), 0)
})

test_that("visped parameter 'compact' works", {
  # Use a specific family known to have > 2 siblings
  # Sire="ZZ5", Dam="Z86" from big_family_size_ped
  target_sire <- "ZZ5"
  target_dam <- "Z86"
  
  # Select the family + parents
  raw_ped <- as.data.table(big_family_size_ped)
  family_inds <- raw_ped[Sire == target_sire & Dam == target_dam, Ind]
  
  # Prepare subset
  subset_ped <- raw_ped[Ind %in% c(target_sire, target_dam, family_inds)]
  
  tidy_fam <- tidyped(subset_ped, addgen = TRUE, addnum = TRUE)
  
  # Without compact
  res_full <- visped(tidy_fam, compact = FALSE, showgraph = FALSE, file = tempfile())
  
  # With compact
  res_compact <- visped(tidy_fam, compact = TRUE, showgraph = FALSE, file = tempfile())
  
  # Check if "compact" nodes exist
  node_types <- igraph::V(res_compact$g)$nodetype
  expect_true("compact" %in% node_types)
  
  # Compact graph should have fewer nodes (real id nodes replaced by one compact node)
  expect_lt(igraph::vcount(res_compact$g), igraph::vcount(res_full$g))
})

test_that("visped parameter 'outline' works", {
  tidy_ped <- tidyped(simple_ped)
  res_outline <- visped(tidy_ped, outline = TRUE, showgraph = FALSE, file = tempfile())
  
  # In outline mode, node sizes are very small
  sizes <- igraph::V(res_outline$g)$size
  # Non-highlighted nodes get size 0.0001
  expect_true(any(sizes <= 0.001))
})

test_that("visped parameter 'cex' works", {
  tidy_ped <- tidyped(simple_ped)
  
  # User specified cex
  my_cex <- 0.75
  res <- visped(tidy_ped, cex = my_cex, showgraph = FALSE, file = tempfile())
  
  # Check label.cex attribute
  real_nodes <- igraph::V(res$g)[nodetype == "real"]
  expect_true(all(real_nodes$label.cex == my_cex))
})

test_that("visped parameter 'showf' displays inbreeding coefficients", {
  # Force add f column if not present or calc it
  tidy_f <- tidyped(simple_ped, inbreed = TRUE)
  
  # Check functionality without conditional dependency if tidyped handles it
  # If tidyped uses nadiv internally, we rely on it being available
  if ("f" %in% colnames(tidy_f)) {
      res <- visped(tidy_f, showf = TRUE, showgraph = FALSE, file = tempfile())
      labels <- igraph::V(res$g)$label
      # Should contain brackets like "[0.003]"
      expect_true(any(grepl("\\[", labels)))
  }
  
  # Test automatic calculation when f is missing
  tidy_no_f <- copy(tidy_f)
  tidy_no_f[, f := NULL]

  expect_message(
    res_auto <- visped(
      tidy_no_f,
      showf = TRUE,
      showgraph = FALSE,
      file = tempfile()
    ),
    "Calculated inbreeding coefficients automatically"
  )

  labels_auto <- igraph::V(res_auto$g)$label
  expect_true(any(grepl("\\[", labels_auto)))
})

test_that("visped parameter 'showf' warns on incomplete pedigrees", {
  tp_full <- tidyped(simple_ped)
  tp_bad <- suppressWarnings(tp_full[Gen > 2])

  expect_false(is_tidyped(tp_bad))
  expect_false("f" %in% names(tp_bad))

  warn_msgs <- character(0)
  res_bad <- withCallingHandlers(
    visped(
      tp_bad,
      showf = TRUE,
      showgraph = FALSE,
      file = tempfile()
    ),
    warning = function(w) {
      warn_msgs <<- c(warn_msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl(
    "cannot be computed automatically because the pedigree is structurally incomplete",
    warn_msgs,
    fixed = TRUE
  )))

  labels_bad <- igraph::V(res_bad$g)$label
  expect_false(any(grepl("\\[", labels_bad)))
})

test_that("visped parameter 'highlight' works with vector", {
  tidy_ped <- tidyped(simple_ped)
  target <- "J5X804"
  
  res <- visped(tidy_ped, highlight = target, showgraph = FALSE, file = tempfile())
  
  hl_nodes <- igraph::V(res$g)[highlighted]
  expect_gt(length(hl_nodes), 0)
  expect_true(target %in% hl_nodes$label)
})

test_that("visped parameter 'highlight' works with list", {
  tidy_ped <- tidyped(simple_ped)
  hl_list <- list(ids = c("J5X804", "J3Y620"), frame.color = "red", color = "yellow")
  
  res <- visped(tidy_ped, highlight = hl_list, showgraph = FALSE, file = tempfile())
  
  hl_nodes <- igraph::V(res$g)[highlighted]
  expect_gt(length(hl_nodes), 0)
})


test_that("visped parameter 'trace' works", {
  tidy_ped <- tidyped(simple_ped)
  target <- "J5X804"
  
  # Trace up
  res_up <- visped(tidy_ped, highlight = target, trace = "up", showgraph = FALSE, file = tempfile())
  
  # Trace down
  res_down <- visped(tidy_ped, highlight = target, trace = "down", showgraph = FALSE, file = tempfile())
  
  # The set of highlighted nodes should differ for this individual (who has parents and offspring)
  h_up <- igraph::V(res_up$g)[highlighted]$label
  h_down <- igraph::V(res_down$g)[highlighted]$label
  
  expect_false(setequal(h_up, h_down))
})

test_that("visped parameter 'file' generates PDF", {
  tidy_ped <- tidyped(simple_ped)
  tmp <- tempfile(fileext = ".pdf")
  
  expect_message(visped(tidy_ped, file = tmp), "Pedigree saved to")
  expect_true(file.exists(tmp))
  unlink(tmp)
})

test_that("visped filters isolated individuals", {
  # Manually construct isolated individual cleanly
  ped_iso <- data.table(
    Ind = c("A", "B", "C", "ISO1"),
    Sire = c(NA_character_, "A", "A", NA_character_),
    Dam = c(NA_character_, NA_character_, NA_character_, NA_character_),
    Sex = c("male", "male", "female", "male")
  )
  
  expect_message(res <- visped(ped_iso, showgraph = FALSE, file = tempfile()), "Removed 1 isolated individuals")
  
  # Check graph does not contain ISO1
  expect_false("ISO1" %in% igraph::V(res$g)$label)
})

test_that("visped parameter 'showgraph' controls plotting", {
    # It's hard to test side effects (plotting) but we can ensure showgraph=FALSE returns silently
    tidy_ped <- tidyped(simple_ped)
    expect_error(visped(tidy_ped, showgraph = FALSE), "Both 'showgraph' and 'file' are disabled")
})

test_that("visped parameter 'pagewidth' works", {
  tidy_ped <- tidyped(simple_ped)
  
  # Default case (should handle standard width)
  res_default <- visped(tidy_ped, showgraph = FALSE, file = tempfile())
  
  # Custom width limit
  # Force a very small width limit to ensure it constrains the output
  small_limit <- 15
  res_small <- visped(tidy_ped, pagewidth = small_limit, showgraph = FALSE, file = tempfile())
  
  expect_lte(res_small$canvas_width, small_limit)
  expect_lte(res_small$canvas_height, small_limit)
  
  # Ensure it doesn't crash on large limits
  res_large <- visped(tidy_ped, pagewidth = 1000, showgraph = FALSE, file = tempfile())
  expect_true(res_large$canvas_width <= 1000)
})

test_that("visped validates new parameters", {
  tidy_ped <- tidyped(simple_ped)
  
  # Validate pagewidth
  expect_error(visped(tidy_ped, pagewidth = -10), "'pagewidth' must be a single positive number")
  expect_error(visped(tidy_ped, pagewidth = "invalid"), "'pagewidth' must be a single positive number")
  
  # Validate symbolsize
  expect_error(visped(tidy_ped, symbolsize = 0), "'symbolsize' must be a single positive number")
  expect_error(visped(tidy_ped, symbolsize = "big"), "'symbolsize' must be a single positive number")
  
  # Validate maxiter
  expect_error(visped(tidy_ped, maxiter = -100), "'maxiter' must be a single positive integer")
  expect_error(visped(tidy_ped, maxiter = "many"), "'maxiter' must be a single positive integer")
  
  # Validate compact
  expect_error(visped(tidy_ped, compact = "yes"), "'compact' must be TRUE or FALSE")
  expect_error(visped(tidy_ped, compact = NA), "'compact' must be TRUE or FALSE")
  expect_error(visped(tidy_ped, compact = NULL), "'compact' must be TRUE or FALSE")
  
  # Validate outline
  expect_error(visped(tidy_ped, outline = 1), "'outline' must be TRUE or FALSE")
  
  # Validate showgraph
  expect_error(visped(tidy_ped, showgraph = "no"), "'showgraph' must be TRUE or FALSE")
  
  # Validate showf
  expect_error(visped(tidy_ped, showf = NA), "'showf' must be TRUE or FALSE")
  
  # Validate cex
  expect_error(visped(tidy_ped, cex = -1), "'cex' must be NULL or a single positive number")
  expect_error(visped(tidy_ped, cex = "small"), "'cex' must be NULL or a single positive number")
  
  # Validate file
  expect_error(visped(tidy_ped, file = TRUE), "'file' must be NULL or a single character string")
  expect_error(visped(tidy_ped, file = NA_character_), "'file' must be NULL or a single character string")

  # Validate highlight
  expect_error(visped(tidy_ped, highlight = 123), "'highlight' must be NULL, a character vector, or a list")

  # Sanitization test
  # Should not crash with NA or empty string
  expect_no_error(suppressMessages(visped(tidy_ped, highlight = c("J5X804", NA, ""), showgraph = FALSE, file = tempfile())))
  expect_no_error(suppressMessages(visped(tidy_ped, highlight = list(ids = c("J5X804", NA, "")), showgraph = FALSE, file = tempfile())))
  
  # Validate trace
  expect_error(visped(tidy_ped, trace = "left"), "'trace' must be TRUE, FALSE, 'up', 'down', or 'all'")
  expect_error(visped(tidy_ped, trace = 1), "'trace' must be TRUE, FALSE, 'up', 'down', or 'all'")
})