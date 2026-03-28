#' Internal layout engine for pedigree visualization
#' @import data.table
#' @importFrom igraph graph_from_data_frame layout_with_sugiyama norm_coords V
#' @importFrom graphics strwidth
#' @keywords internal
prepare_ped_graph <- function(ped, compact = FALSE, outline = FALSE, cex = NULL, 
                              highlight = NULL, trace = FALSE, showf = FALSE, pagewidth = 200, symbolsize = 1, maxiter = 1000, ...) {
  ped_new <- copy(ped)

  # Check and tidyped if necessary
  if (!"Gen" %in% colnames(ped_new) || !all(c("IndNum", "SireNum", "DamNum") %in% colnames(ped_new))) {
    ped_new <- tidyped(ped_new, addgen = TRUE, addnum = TRUE)
  }
  
  # Remove isolated individuals early to optimize performance and simplify graph conversion
  if ("Gen" %in% colnames(ped_new)) {
    n_iso <- sum(ped_new$Gen == 0, na.rm = TRUE)
    if (n_iso > 0) {
      message(sprintf("Note: Removed %d isolated individuals (no parents, no progeny) from the plot.", n_iso))
      ped_new <- ped_new[Gen != 0]
    }
  }

  # Reserved digits
  fixed_digits <- 7
  old_digits <- getOption("digits")
  options(digits = 20)
  on.exit(options(digits = old_digits), add = TRUE)

  ped_igraph_data <- ped2igraph(ped_new, compact, highlight, trace, showf)
  ped_igraph <- ped_igraph_data
  real_node <- ped_igraph$node[nodetype %in% c("real", "compact")]
  
  if (nrow(real_node) == 0) {
    stop("Pedigree contains no individuals to plot.")
  }

  gen_node_num <- real_node[, .N, by = gen]
  gen_max_size <- max(gen_node_num$N, na.rm = TRUE)

  pdf_max_width_soft <- pagewidth
  s_size <- symbolsize
  cexs <- seq(from = 0.1, to = 1, by = 0.05)
  best_cex <- 0
  max_strwidth_label <- real_node[which.max(strwidth(real_node$label, cex = 1, units = "inches")), label]
  
  for (i in length(cexs):1) {
    label_max_width <- max(strwidth(max_strwidth_label, cex = cexs[i], units = "inches"), na.rm = TRUE)
    if (gen_max_size <= 16 & label_max_width < 0.8) {
      label_max_width <- 0.8
    }
    if ((label_max_width * s_size * gen_max_size) < pdf_max_width_soft) {
      best_cex <- cexs[i] * 0.65
      break
    }
  }

  if (!outline & best_cex == 0 & is.null(cex)) {
    stop(sprintf("Too many individuals per generation; cannot fit within %.1f inches. Use compact=TRUE, outline=TRUE, or increase pagewidth.", pdf_max_width_soft))
  }

  if (!outline && !is.null(cex)) {
    label_max_width <- max(strwidth(max_strwidth_label, cex = cex, units = "inches"), na.rm = TRUE)
    if (gen_max_size <= 16 && label_max_width < 0.8) {
      label_max_width <- 0.8
    }
    if ((label_max_width * s_size * gen_max_size) >= pdf_max_width_soft) {
      warning(sprintf("Provided cex likely exceeds PDF width (%.1f inches); labels may clip. Try compact=TRUE, outline=TRUE, or increase pagewidth.", pdf_max_width_soft), call. = FALSE)
    }
  }

  hgap <- round(1 / gen_max_size, 8)
  gen_num <- max(real_node$gen, na.rm = TRUE)
  max_layer <- max(ped_igraph$node$layer, na.rm = TRUE)
  
  # Create a minimal graph for layout calculation to save memory and time
  g_layout <- graph_from_data_frame(ped_igraph$edge[, .(from, to)], directed = TRUE, vertices = ped_igraph$node[, .(id)])
  # match names to get layers - V(g_layout)$name exists by default from graph_from_data_frame
  layer_idx <- ped_igraph$node$layer[match(v_names <- V(g_layout)$name, as.character(ped_igraph$node$id))]
  layer_idx <- max_layer - layer_idx + 1
  
  l <- layout_with_sugiyama(g_layout, layers = layer_idx, hgap = hgap, maxiter = maxiter, attributes = "all")$layout
  l <- norm_coords(l, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
  
  # Map coordinates back to ped_igraph$node based on names
  # The layout l corresponds to V(g_layout)
  layout_dt <- data.table(id = as.integer(v_names), x = l[, 1], y = l[, 2])
  ped_igraph$node[layout_dt, on = "id", `:=`(x = i.x, y = i.y)]

  real_node <- ped_igraph$node[nodetype %in% c("real", "compact")]
  for (i in 1:gen_num) {
    v_rank <- rank(real_node[gen == i, x], na.last = TRUE, ties.method = "first")
    x_sorted <- round(sort(real_node[gen == i, x]), fixed_digits)
    x_new <- repeloverlap(x_sorted)
    real_node[gen == i, x := x_new[v_rank]]
  }
  ped_igraph$node[nodetype %in% c("real", "compact")] <- real_node

  if (gen_max_size >= 2) {
    min_gap <- 0.8 / gen_max_size
    fam_gap <- min_gap * 1.2
    
    for (i in seq_len(gen_num)) {
      current_nodes <- real_node[gen == i]
      if (nrow(current_nodes) == 0) next
      
      if (any(is.na(current_nodes$familylabel))) {
        current_nodes[is.na(familylabel), familylabel := paste0("u_", id)]
      }
      
      parent_coords <- real_node[, .(Ind, px = x)]
      current_nodes[parent_coords, on = .(sirelabel = Ind), sire_x := i.px]
      current_nodes[parent_coords, on = .(damlabel = Ind), dam_x := i.px]
      current_nodes[, parent_mid := rowMeans(.SD, na.rm = TRUE), .SDcols = c("sire_x", "dam_x")]
      current_nodes[is.na(parent_mid), parent_mid := x]
      
      fam_stats <- current_nodes[, .(
        N = .N,
        target_x = mean(parent_mid, na.rm=TRUE),
        original_x_mean = mean(x, na.rm=TRUE)
      ), by = familylabel]
      
      fam_stats[, half_width := ((N - 1) * min_gap) / 2]
      setorder(fam_stats, target_x, original_x_mean)
      
      n_fams <- nrow(fam_stats)
      if (n_fams > 0) {
        fam_stats[, fwd_pos := target_x]
        fam_stats[, bwd_pos := target_x]
        
        if (n_fams > 1) {
          for (k in 2:n_fams) {
            req_dist <- fam_stats$half_width[k-1] + fam_gap + fam_stats$half_width[k]
            if (fam_stats$fwd_pos[k] < fam_stats$fwd_pos[k-1] + req_dist) {
              fam_stats$fwd_pos[k] <- fam_stats$fwd_pos[k-1] + req_dist
            }
          }
          for (k in (n_fams-1):1) {
            req_dist <- fam_stats$half_width[k] + fam_gap + fam_stats$half_width[k+1]
            if (fam_stats$bwd_pos[k] > fam_stats$bwd_pos[k+1] - req_dist) {
              fam_stats$bwd_pos[k] <- fam_stats$bwd_pos[k+1] - req_dist
            }
          }
        }
        
        fam_stats[, final_center := (fwd_pos + bwd_pos) / 2]
        current_nodes[fam_stats, on = "familylabel", fam_center := i.final_center]
        setorder(current_nodes, familylabel, x) 
        current_nodes[, sib_offset := (seq_len(.N) - 1 - (.N - 1)/2) * min_gap, by = familylabel]
        real_node[current_nodes, on = "id", x := i.fam_center + i.sib_offset]
      }
    }
    
    rx <- range(real_node$x, na.rm = TRUE)
    if (diff(rx) > 1e-9) {
      real_node[, x := (x - rx[1]) / diff(rx)]
    } else {
      real_node[, x := 0.5]
    }
    ped_igraph$node[nodetype %in% c("real", "compact")] <- real_node
  }

  virtual_node <- ped_igraph$node[nodetype %in% c("virtual")]
  if (nrow(virtual_node) > 0) {
    real_family_min_x <- real_node[, .(minx = min(x, na.rm = TRUE)), by = .(gen, familylabel)]
    virtual_node[real_family_min_x, x := i.minx, on = .(gen, familylabel)]
  }
  ped_igraph$node[nodetype %in% c("virtual")] <- virtual_node

  # Update layout matrix with adjusted x coordinates for correct canvas width calculation
  l[, 1] <- ped_igraph$node[, x]

  s_size <- if (!is.null(symbolsize)) symbolsize else 1
  calc_node_width <- label_max_width * s_size
  node_width_s <- label_max_width * s_size
  
  if (!outline) {
    if (gen_max_size >= 2) {
      x_stats_gen <- real_node[, .(.N, range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE)), by = gen]
      x_stats_gen[range > 0 & N > 1, ":="(meanspace = range / (N - 1))]
      min_node_space <- min(x_stats_gen[meanspace > 0, meanspace], na.rm = TRUE)
      f <- if (max(x_stats_gen$N) <= 16) 3 * round(node_width_s / min_node_space, 8) else round(node_width_s / min_node_space, 8)
    } else {
      f <- 1
    }
    canvas_width_s <- max(f * l[, 1], na.rm = TRUE) - min(f * l[, 1], na.rm = TRUE) + 6 * node_width_s
  } else {
    canvas_width_s <- min(label_max_width * gen_max_size, pdf_max_width_soft)
    node_width_s <- 0.0001
    if (0.0001 * gen_max_size > pdf_max_width_soft) stop("Outline view unavailable.")
  }

  canvas_width_s <- max(10, min(pdf_max_width_soft, canvas_width_s))
  min_layer_height <- if (outline) 0.1 else 0.8
  effective_layer_height <- max(calc_node_width, min_layer_height)
  canvas_height <- gen_num * effective_layer_height + 3 * effective_layer_height
  canvas_height <- max(8, canvas_height, canvas_width_s * 0.618)

  # Prepare final attributes in data.table BEFORE creating the graph
  ped_igraph$node[, `:=`(size = 0.001, label.cex = 1.0)]
  
  node_size_val <- round(node_width_s * 100 / canvas_width_s, 8)
  current_cex <- if (is.null(cex)) best_cex else cex
  
  if (outline) {
    # Outline mode: icons are tiny, no labels except highlighted
    ped_igraph$node[nodetype %in% c("real", "compact"), size := 0.0001]
    is_h <- !is.na(ped_igraph$node$highlighted) & ped_igraph$node$highlighted
    ped_igraph$node[!is_h, label := ""]
    
    if (any(is_h)) {
      highlight_node_size <- round(label_max_width * 100 / canvas_width_s, 8)
      ped_igraph$node[is_h & nodetype %in% c("real", "compact"), `:=`(
        size = highlight_node_size,
        label.cex = if (is.null(cex)) (if (best_cex > 0) best_cex else 0.6) else cex
      )]
    }
  } else {
    ped_igraph$node[nodetype %in% c("real", "compact"), `:=`(
      size = node_size_val,
      label.cex = current_cex
    )]
  }

  # Compute generation info for labels (using final layout coordinates)
  gen_info <- ped_igraph$node[nodetype %in% c("real", "compact"),
                              .(y = mean(y, na.rm = TRUE), N = .N), by = gen]
  setorder(gen_info, gen)
  
  # Final graph creation with all attributes
  # Ensure the order of columns in edge and node are correct for igraph
  g <- graph_from_data_frame(ped_igraph$edge, directed = TRUE, ped_igraph$node)
  max_node_size <- max(ped_igraph$node$size, na.rm = TRUE)

  return(list(
    g = g,
    layout = as.matrix(ped_igraph$node[, .(x, y)]),
    canvas_width = canvas_width_s,
    canvas_height = canvas_height,
    node_size = max_node_size,
    best_cex = best_cex,
    gen_info = gen_info
  ))
}

#' Repel overlapping nodes on the x-axis
#' @param x A numeric vector of x positions.
#' @return A numeric vector of unique x positions.
#' @keywords internal
repeloverlap <- function(x) {
  n <- length(x)
  if (n <= 1 || anyDuplicated(x) == 0) return(x)
  
  x_dt <- as.data.table(x)
  x_counts <- x_dt[, .N, by = x]
  dup_info <- x_counts[N > 1]
  if (nrow(dup_info) == 0 || nrow(x_counts) == 1) return(x)
  
  setorder(dup_info, x)
  unique_pos <- sort(x_counts[, x])
  n_unique <- length(unique_pos)
  n_dup <- nrow(dup_info)
  
  max_new_vals <- sum(dup_info$N) - n_dup
  result <- numeric(n_unique + max_new_vals)
  result[seq_len(n_unique)] <- unique_pos
  result_idx <- n_unique
  
  for (i in seq_len(n_dup)) {
    dup_val <- dup_info$x[i]
    n_copies <- dup_info$N[i]
    idx <- match(dup_val, unique_pos)
    
    if (i == n_dup - 1 && n_dup >= 2) {
      next_dup_val <- dup_info$x[i + 1]
      next_idx <- match(next_dup_val, unique_pos)
      if (next_idx == idx + 1) {
        next_n_copies <- dup_info$N[i + 1]
        total_copies <- n_copies + next_n_copies - 1
        gap <- (next_dup_val - dup_val) / total_copies
        n_new1 <- n_copies - 1
        n_new2 <- next_n_copies - 1
        result[(result_idx + 1):(result_idx + n_new1)] <- dup_val + gap * seq_len(n_new1)
        result[(result_idx + n_new1 + 1):(result_idx + n_new1 + n_new2)] <- 
          dup_val + gap * seq(n_copies, total_copies - 1)
        result_idx <- result_idx + n_new1 + n_new2
        break
      }
    }
    
    n_new <- n_copies - 1
    if (n_new > 0) {
      if (idx < n_unique) {
        next_pos <- unique_pos[idx + 1]
        gap <- (next_pos - dup_val) / n_copies
        result[(result_idx + 1):(result_idx + n_new)] <- dup_val + gap * seq_len(n_new)
      } else if (idx > 1) {
        prev_pos <- unique_pos[idx - 1]
        gap <- (dup_val - prev_pos) / n_copies
        result[(result_idx + 1):(result_idx + n_new)] <- prev_pos + gap * seq_len(n_new)
      }
      result_idx <- result_idx + n_new
    }
  }
  return(result[seq_len(result_idx)])
}
