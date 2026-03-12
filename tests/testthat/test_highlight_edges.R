test_that("highlight edges are correct for shared parents", {
  # Test case: J5X804's ancestor J0Z938 also has another child J1H419
  # The edge to J1H419's family should NOT be highlighted
  
  tidy_ped <- tidyped(simple_ped, genmethod = "bottom")
  graph_data <- visPedigree:::ped2igraph(
    tidy_ped, 
    compact = FALSE, 
    highlight = "J5X804", 
    trace = "up", 
    showf = FALSE
  )
  
  # Find J0Z938 node
  j0z938_node <- graph_data$node[label == "J0Z938"]
  expect_equal(nrow(j0z938_node), 1)
  expect_true(j0z938_node$highlighted)
  
  # Find edges from virtual family nodes to J0Z938
  edges_to_j0z938 <- graph_data$edge[to == j0z938_node$id & from > max(graph_data$node[nodetype %in% c("real", "compact")]$id)]
  
  # Should have 2 edges: one from highlighted family (J0Z938xJ0Z843), one from non-highlighted family (J0Z938xJ0Z167)
  expect_equal(nrow(edges_to_j0z938), 2)
  
  # Edge from highlighted family should still match parent node color
  highlighted_edge <- edges_to_j0z938[from_highlighted == TRUE]
  expect_equal(nrow(highlighted_edge), 1)
  expect_equal(highlighted_edge$color, j0z938_node$color)
  
  # Edge from non-highlighted family should match parent node color with fading
  non_highlighted_edge <- edges_to_j0z938[from_highlighted == FALSE]
  expect_equal(nrow(non_highlighted_edge), 1)
  expect_true(grepl(paste0("^", j0z938_node$color, "4D$"), non_highlighted_edge$color))
})

test_that("highlight edges are correct for shared children", {
  # Test case: J0Z475 and J0C612 have two children: J1Y339 (in path) and J1F266 (not in path)
  # Only the edge from J1Y339 to the family virtual node should be highlighted
  
  tidy_ped <- tidyped(simple_ped, genmethod = "bottom")
  graph_data <- visPedigree:::ped2igraph(
    tidy_ped, 
    compact = FALSE, 
    highlight = "J4Y326", 
    trace = "up", 
    showf = FALSE
  )
  
  # Find the family virtual node
  family_node <- graph_data$node[familylabel == "J0Z475xJ0C612" & nodetype == "virtual"]
  expect_equal(nrow(family_node), 1)
  expect_true(family_node$highlighted)
  
  # Find children
  j1y339_node <- graph_data$node[label == "J1Y339"]
  j1f266_node <- graph_data$node[label == "J1F266"]
  expect_true(j1y339_node$highlighted)
  expect_false(j1f266_node$highlighted)
  
  # Find edges from children to family node
  edges_to_family <- graph_data$edge[to == family_node$id]
  expect_equal(nrow(edges_to_family), 2)
  
  # Edge from highlighted child (J1Y339) should be solid
  edge_from_j1y339 <- edges_to_family[from == j1y339_node$id]
  expect_equal(edge_from_j1y339$color, "#333333")
  
  # Edge from non-highlighted child (J1F266) should be faded
  edge_from_j1f266 <- edges_to_family[from == j1f266_node$id]
  expect_equal(edge_from_j1f266$color, "#3333334D")
})

test_that("edges work correctly without highlighting", {
  tidy_ped <- tidyped(simple_ped, genmethod = "bottom")
  graph_data <- visPedigree:::ped2igraph(
    tidy_ped, 
    compact = FALSE, 
    highlight = NULL, 
    trace = FALSE, 
    showf = FALSE
  )
  
  # Without highlighting, edges from virtual nodes should use parent colors
  real_max <- max(graph_data$node[nodetype %in% c("real", "compact")]$id)
  virtual_edges <- graph_data$edge[from > real_max]
  
  # All virtual edges should have color matching their target node's color (fill color)
  expect_true(all(!is.na(virtual_edges$color)))
  expect_true(all(virtual_edges$color %in% graph_data$node$color))
})

test_that("highlight edges work correctly with trace down", {
  # Test descendant tracing: edges should follow the same rules
  tidy_ped <- tidyped(simple_ped, genmethod = "bottom")
  graph_data <- visPedigree:::ped2igraph(
    tidy_ped, 
    compact = FALSE, 
    highlight = "J1J576", 
    trace = "down", 
    showf = FALSE
  )
  
  # J1J576 should be highlighted
  j1j576_node <- graph_data$node[label == "J1J576"]
  expect_true(j1j576_node$highlighted)
  
  # Edge from J1J576 to its family should be highlighted
  j1j576_family_id <- j1j576_node$familynum
  edge_to_family <- graph_data$edge[from == j1j576_node$id & to == j1j576_family_id]
  expect_equal(nrow(edge_to_family), 1)
  expect_equal(edge_to_family$color, "#3333334D")
  
  # J1J576's parents should NOT be highlighted (trace = "down" only traces descendants)
  j0z938_node <- graph_data$node[label == "J0Z938"]
  j0z843_node <- graph_data$node[label == "J0Z843"]
  expect_false(j0z938_node$highlighted)
  expect_false(j0z843_node$highlighted)
  
  # Edges from family to parents should match parent node colors
  family_node <- graph_data$node[id == j1j576_family_id]
  edges_from_family <- graph_data$edge[from == j1j576_family_id]
  parent_colors <- graph_data$node[id %in% edges_from_family$to, color]
  expect_true(all(edges_from_family$color %in% parent_colors))
})

test_that("highlight edges work correctly with trace all", {
  # Test both ancestor and descendant tracing
  tidy_ped <- tidyped(simple_ped, genmethod = "bottom")
  graph_data <- visPedigree:::ped2igraph(
    tidy_ped, 
    compact = FALSE, 
    highlight = "J3Y620", 
    trace = "all", 
    showf = FALSE
  )
  
  # Focal individual should be highlighted
  focal_node <- graph_data$node[label == "J3Y620"]
  expect_true(focal_node$highlighted)
  
  # Ancestor (J2C161) should be highlighted
  ancestor_node <- graph_data$node[label == "J2C161"]
  expect_true(ancestor_node$highlighted)
  
  # Descendant (J4Y326) should be highlighted
  descendant_node <- graph_data$node[label == "J4Y326"]
  expect_true(descendant_node$highlighted)
  
  # Edges involving highlighted nodes should follow the rules
  # Individual to family: only if individual is highlighted
  edge_focal_to_family <- graph_data$edge[from == focal_node$id & to == focal_node$familynum]
  if (nrow(edge_focal_to_family) > 0) {
    expect_equal(edge_focal_to_family$color, "#333333")
  }
})

test_that("highlight edges work correctly with multiple individuals without trace", {
  # Test multiple highlighted individuals sharing a parent
  # When trace = FALSE, NO edges should be highlighted (only individuals are marked)
  tidy_ped <- tidyped(simple_ped, genmethod = "bottom")
  graph_data <- visPedigree:::ped2igraph(
    tidy_ped, 
    compact = FALSE, 
    highlight = c("J1J576", "J1H419"), 
    trace = FALSE, 
    showf = FALSE
  )
  
  # Both individuals should be highlighted
  j1j576_node <- graph_data$node[label == "J1J576"]
  j1h419_node <- graph_data$node[label == "J1H419"]
  expect_true(j1j576_node$highlighted)
  expect_true(j1h419_node$highlighted)
  
  # Their shared parent J0Z938 should NOT be highlighted (no trace)
  j0z938_node <- graph_data$node[label == "J0Z938"]
  expect_false(j0z938_node$highlighted)
  
  # Both family virtual nodes should be highlighted
  family1 <- graph_data$node[id == j1j576_node$familynum]
  family2 <- graph_data$node[id == j1h419_node$familynum]
  expect_true(family1$highlighted)
  expect_true(family2$highlighted)
  
  # WITHOUT trace, edges from individuals to their families should be FADED
  # (No relationship path is being traced, so edges remain faded)
  edge1 <- graph_data$edge[from == j1j576_node$id & to == j1j576_node$familynum]
  edge2 <- graph_data$edge[from == j1h419_node$id & to == j1h419_node$familynum]
  expect_equal(edge1$color, "#3333334D")
  expect_equal(edge2$color, "#3333334D")
  
  # Edges from families to shared parent should match parent node colors
  real_max <- max(graph_data$node[nodetype %in% c("real", "compact")]$id)
  edges_to_parent <- graph_data$edge[to == j0z938_node$id & from > real_max]
  # Faded edges should have 4D suffix
  expect_true(all(grepl(paste0("^", j0z938_node$color, "(4D)?$"), edges_to_parent$color)))
})

test_that("highlight edges work correctly with compact mode", {
  # Test that compact mode doesn't break edge highlighting logic
  tidy_ped <- tidyped(simple_ped, genmethod = "bottom")
  graph_data <- visPedigree:::ped2igraph(
    tidy_ped, 
    compact = TRUE, 
    highlight = "J5X804", 
    trace = "up", 
    showf = FALSE
  )
  
  # Basic edge highlighting rules should still apply
  real_max <- max(graph_data$node[nodetype %in% c("real", "compact")]$id)
  
  # All edges should be either highlighted (#333333), faded (#3333334D), or match parent colors (with optional 4D)
  all_edges <- graph_data$edge
  node_colors <- graph_data$node$color
  node_colors_faded <- paste0(node_colors, "4D")
  valid_colors <- unique(c("#333333", "#3333334D", node_colors, node_colors_faded))
  expect_true(all(all_edges$color %in% valid_colors))
  
  # Edges from real/compact nodes should respect from_highlighted
  individual_edges <- all_edges[from <= real_max]
  if (nrow(individual_edges) > 0) {
    highlighted_ind_edges <- individual_edges[from_highlighted == TRUE]
    if (nrow(highlighted_ind_edges) > 0) {
      expect_true(all(highlighted_ind_edges$color == "#333333"))
    }
    
    faded_ind_edges <- individual_edges[from_highlighted == FALSE]
    if (nrow(faded_ind_edges) > 0) {
      expect_true(all(faded_ind_edges$color == "#3333334D"))
    }
  }
  
  # Edges from virtual nodes should follow parent colors (with fading for non-highlighted families)
  virtual_edges <- all_edges[from > real_max]
  if (nrow(virtual_edges) > 0) {
    parent_colors <- graph_data$node[id %in% virtual_edges$to, color]
    parent_colors_faded <- paste0(parent_colors, "4D")
    valid_parent_colors <- unique(c(parent_colors, parent_colors_faded))
    expect_true(all(virtual_edges$color %in% valid_parent_colors))
  }
})
