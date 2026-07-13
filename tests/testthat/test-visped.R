library(testthat)
library(data.table)
library(visPedigree)

make_sex_shape_ped <- function() {
  tidyped(
    data.table(
      Ind = c("M", "F", "H", "U", "HC"),
      Sire = c(NA, NA, NA, "M", "H"),
      Dam = c(NA, NA, NA, "F", "H")
    ),
    selfing = TRUE
  )
}

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

test_that("visped writes SVG output when requested", {
  tidy_ped <- tidyped(simple_ped, addgen = TRUE, addnum = TRUE)
  svg_file <- tempfile(fileext = ".svg")

  res <- visped(tidy_ped, showgraph = FALSE, file = svg_file)

  expect_s3_class(res$g, "igraph")
  expect_true(file.exists(svg_file))
  expect_gt(file.info(svg_file)$size, 0)
})

test_that("visped supports custom generation labels in displayed-layer order", {
  tidy_ped <- tidyped(simple_ped, addgen = TRUE, addnum = TRUE)
  default_info <- visPedigree:::prepare_visped_geninfo(
    data.table(gen = c(1L, 3L, 5L), y = c(0.1, 0.5, 0.9))
  )
  custom_info <- visPedigree:::prepare_visped_geninfo(
    data.table(gen = c(1L, 3L, 5L), y = c(0.1, 0.5, 0.9)),
    c("Founders", "Cycle 1", "Current")
  )

  expect_equal(default_info$label, c("G1", "G3", "G5"))
  expect_equal(custom_info$label, c("Founders", "Cycle 1", "Current"))

  custom_labels <- paste0("Layer ", seq_len(length(unique(tidy_ped$Gen))))
  expect_silent(
    suppressMessages(visped(
      tidy_ped,
      genlab = custom_labels,
      showgraph = FALSE,
      file = tempfile(fileext = ".pdf")
    ))
  )
})

test_that("visped validates custom generation labels", {
  tidy_ped <- tidyped(simple_ped, addgen = TRUE, addnum = TRUE)

  expect_error(
    visped(tidy_ped, genlab = c("one", "two"), showgraph = FALSE, file = tempfile()),
    "one label for each displayed generation"
  )
  expect_error(
    visped(tidy_ped, genlab = c("1" = "Founders"), showgraph = FALSE, file = tempfile()),
    "unnamed character vector"
  )
  expect_error(
    visped(tidy_ped, genlab = character(), showgraph = FALSE, file = tempfile()),
    "TRUE, FALSE, or a non-empty unnamed character vector"
  )
  expect_error(
    visped(tidy_ped, genlab = 1, showgraph = FALSE, file = tempfile()),
    "TRUE, FALSE, or a non-empty unnamed character vector"
  )
})

test_that("visped reports the effective generation-label cex", {
  tidy_ped <- tidyped(simple_ped, addgen = TRUE, addnum = TRUE)

  expect_message(
    visped(
      tidy_ped,
      genlab = TRUE,
      genlabcex = 1.2,
      showgraph = FALSE,
      file = tempfile(fileext = ".pdf")
    ),
    "Generation label cex: 1.2"
  )
})

test_that("visped supports custom individual labels", {
  tidy_ped <- tidyped(simple_ped, addgen = TRUE, addnum = TRUE)
  tidy_ped[, DisplayID := paste0("node_", seq_len(.N))]

  res <- visped(tidy_ped, labelvar = "DisplayID", showgraph = FALSE, file = tempfile())
  vertices <- as.data.table(igraph::as_data_frame(res$g, what = "vertices"))
  real_nodes <- vertices[nodetype == "real"]

  expected <- tidy_ped$DisplayID[match(real_nodes$Ind, tidy_ped$Ind)]
  expect_equal(real_nodes$label, expected)
})

test_that("visped accepts row-aligned labels without losing ID alignment", {
  raw_ped <- copy(as.data.table(simple_ped))
  raw_ped <- raw_ped[rev(seq_len(nrow(raw_ped)))]
  original_names <- names(raw_ped)
  raw_ids <- as.character(raw_ped[[1L]])
  display_labels <- paste0("node_", raw_ids)
  display_labels[1:2] <- c(NA_character_, "")
  focal <- "J5X804"

  res <- visped(
    raw_ped,
    labelvar = display_labels,
    highlight = focal,
    trace = "all",
    showgraph = FALSE,
    file = tempfile()
  )
  vertices <- as.data.table(igraph::as_data_frame(res$g, what = "vertices"))
  real_nodes <- vertices[nodetype == "real"]
  expected <- display_labels[match(real_nodes$Ind, raw_ids)]
  expected[is.na(expected) | expected == ""] <-
    real_nodes$Ind[is.na(expected) | expected == ""]

  expect_equal(real_nodes$label, expected)
  expect_true(real_nodes[Ind == focal, highlighted])
  expect_identical(names(raw_ped), original_names)
})

test_that("visped rejects missing custom label columns", {
  tidy_ped <- tidyped(simple_ped, addgen = TRUE, addnum = TRUE)

  expect_error(
    visped(tidy_ped, labelvar = "MissingLabel", showgraph = FALSE, file = tempfile()),
    "specified by 'labelvar' was not found"
  )
  expect_error(
    visped(tidy_ped, labelvar = c("one", "two")),
    "one label per pedigree row"
  )
  expect_error(
    visped(tidy_ped, labelvar = seq_len(nrow(tidy_ped))),
    "must be NULL"
  )
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
  res_compact <- visped(
    tidy_fam,
    compact = TRUE,
    labelvar = paste0("node_", tidy_fam$Ind),
    showgraph = FALSE,
    file = tempfile()
  )
  
  # Check if "compact" nodes exist
  node_types <- igraph::V(res_compact$g)$nodetype
  expect_true("compact" %in% node_types)
  compact_labels <- igraph::V(res_compact$g)[nodetype == "compact"]$label
  real_labels <- igraph::V(res_compact$g)[nodetype == "real"]$label
  expect_true(all(grepl("^FS\u00d7[0-9]+$", compact_labels)))
  expect_true(all(grepl("^node_", real_labels)))
  expect_true(all(igraph::V(res_compact$g)[nodetype == "compact"]$shape == "rectangle"))
  expect_equal(
    igraph::V(res_compact$g)[nodetype == "compact"]$size2,
    igraph::V(res_compact$g)[nodetype == "compact"]$size
  )
  expect_true(all(
    igraph::V(res_compact$g)[nodetype == "compact"]$color == "#9cb383"
  ))
  expect_true(all(
    igraph::V(res_compact$g)[nodetype == "compact"]$frame.color == "#5f7650"
  ))
  
  # Compact graph should have fewer nodes (real id nodes replaced by one compact node)
  expect_lt(igraph::vcount(res_compact$g), igraph::vcount(res_full$g))
})

test_that("visped shapeby sex is the default and preserves sex colors", {
  tidy_ped <- make_sex_shape_ped()
  res_default <- visped(tidy_ped, showgraph = FALSE, file = tempfile())
  res_sex <- visped(
    tidy_ped,
    shapeby = "sex",
    showgraph = FALSE,
    file = tempfile()
  )

  default_nodes <- as.data.table(igraph::as_data_frame(res_default$g, what = "vertices"))
  sex_nodes <- as.data.table(igraph::as_data_frame(res_sex$g, what = "vertices"))
  expect_equal(default_nodes$shape, sex_nodes$shape)
  expect_equal(sex_nodes[Ind == "F", shape], "circle")
  expect_equal(sex_nodes[Ind == "M", shape], "square")
  expect_equal(sex_nodes[Ind == "U", shape], "visped_diamond")
  expect_equal(sex_nodes[Ind == "H", shape], "visped_hexagon")
  expect_equal(sex_nodes[Ind == "M", color], "#119ecc")
  expect_equal(sex_nodes[Ind == "F", color], "#f4b131")
  expect_equal(sex_nodes[Ind == "H", color], "#26a69a")
  expect_equal(sex_nodes[Ind == "U", color], "#d9d9d9")
  expect_equal(sex_nodes[Ind == "U", frame.color], "#777777")
})

test_that("visped shapeby role keeps individual records as circles", {
  tidy_ped <- make_sex_shape_ped()
  res <- visped(
    tidy_ped,
    shapeby = "role",
    showgraph = FALSE,
    file = tempfile()
  )
  nodes <- as.data.table(igraph::as_data_frame(res$g, what = "vertices"))

  expect_true(all(nodes[nodetype == "real", shape] == "circle"))
  expect_equal(nodes[Ind == "M", color], "#119ecc")
  expect_equal(nodes[Ind == "F", color], "#f4b131")
  expect_equal(nodes[Ind == "H", color], "#26a69a")
  expect_equal(nodes[Ind == "U", color], "#d9d9d9")
  expect_equal(nodes[Ind == "U", frame.color], "#777777")
})

test_that("visped shapeby sex maps all supported sex classes", {
  tidy_ped <- make_sex_shape_ped()
  res <- visped(
    tidy_ped,
    shapeby = "sex",
    showgraph = FALSE,
    file = tempfile()
  )
  nodes <- as.data.table(igraph::as_data_frame(res$g, what = "vertices"))

  expect_equal(nodes[Ind == "F", shape], "circle")
  expect_equal(nodes[Ind == "M", shape], "square")
  expect_equal(nodes[Ind == "U", shape], "visped_diamond")
  expect_equal(nodes[Ind == "H", shape], "visped_hexagon")
  expect_true(all(c("visped_diamond", "visped_hexagon") %in% igraph::shapes()))
  expect_silent(visPedigree:::register_visped_shapes(c("diamond", "hexagon")))
  expect_error(
    visped(tidy_ped, shapeby = "unsupported"),
    "'arg' should be one of"
  )
})

test_that("custom shapes skip zero-size edge-only rendering pass", {
  zero_params <- function(type, name) {
    switch(
      name,
      color = "#000000",
      frame.color = "#000000",
      frame.width = 0,
      size = 0
    )
  }
  coords <- matrix(c(0.5, 0.5), nrow = 1L)

  expect_silent(visPedigree:::make_visped_polygon(4L)(coords, params = zero_params))
  expect_silent(visPedigree:::make_visped_polygon(6L)(coords, params = zero_params))
})

test_that("compact labels cannot be confused with numeric individual IDs", {
  sibling_ids <- paste0("S", seq_len(10L))
  raw_ped <- data.table(
    Ind = c("10", "F", "P1", "P2", sibling_ids, "C"),
    Sire = c(NA, NA, NA, NA, rep("P1", 10L), "10"),
    Dam = c(NA, NA, NA, NA, rep("P2", 10L), "F")
  )
  res <- visped(
    tidyped(raw_ped),
    compact = TRUE,
    showgraph = FALSE,
    file = tempfile()
  )
  nodes <- as.data.table(igraph::as_data_frame(res$g, what = "vertices"))

  expect_true("10" %in% nodes[nodetype == "real", label])
  expect_true("FS\u00d710" %in% nodes[nodetype == "compact", label])
  expect_equal(nodes[label == "FS\u00d710", shape], "rectangle")
})

test_that("highlighted and parent individuals are not compacted", {
  raw_ped <- data.table(
    Ind = c("P1", "P2", "A", "B", "C", "D"),
    Sire = c(NA, NA, "P1", "P1", "P1", "A"),
    Dam = c(NA, NA, "P2", "P2", "P2", "B")
  )
  tidy_ped <- tidyped(raw_ped)
  res <- visped(
    tidy_ped,
    compact = TRUE,
    highlight = "C",
    showgraph = FALSE,
    file = tempfile()
  )
  nodes <- as.data.table(igraph::as_data_frame(res$g, what = "vertices"))

  expect_true(all(c("A", "B", "C") %in% nodes[nodetype == "real", Ind]))
  expect_false(any(nodes[nodetype == "compact", label] == "FS\u00d73"))
})

test_that("sex-based custom shapes render to PDF and SVG", {
  tidy_ped <- make_sex_shape_ped()
  for (extension in c(".pdf", ".svg")) {
    output <- tempfile(fileext = extension)
    expect_silent(
      suppressMessages(visped(
        tidy_ped,
        shapeby = "sex",
        showgraph = FALSE,
        file = output
      ))
    )
    expect_true(file.exists(output))
    expect_gt(file.info(output)$size, 0)
  }
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

  # Validate genlabcex
  expect_error(visped(tidy_ped, genlabcex = 0), "'genlabcex' must be NULL or a single positive number")
  expect_error(visped(tidy_ped, genlabcex = "large"), "'genlabcex' must be NULL or a single positive number")
})

test_that("visped parameter 'genlabcex' works independently of label cex", {
  tidy_ped <- tidyped(simple_ped)
  tmp <- tempfile(fileext = ".pdf")

  expect_no_error(
    visped(
      tidy_ped,
      showgraph = FALSE,
      file = tmp,
      genlab = TRUE,
      genlabcex = 1.4,
      cex = 0.4
    )
  )

  expect_true(file.exists(tmp))
  unlink(tmp)
})
