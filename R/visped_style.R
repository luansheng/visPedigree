#' Fade colors by appending a reduced alpha value
#'
#' Converts any R color specification to `#RRGGBB4D` form.
#' Handles hex colors (`#RRGGBB`, `#RRGGBBAA`) and named colors (e.g. `"red"`).
#'
#' @param x Character vector of colors.
#' @return Character vector of faded hex colors.
#' @keywords internal
fade_cols <- function(x) {
  vapply(x, function(col) {
    nc <- nchar(col)
    if (nc == 7 && substring(col, 1, 1) == "#") {
      # #RRGGBB -> append 4D
      paste0(col, "4D")
    } else if (nc == 9 && substring(col, 1, 1) == "#") {
      # #RRGGBBAA -> replace AA with 4D
      paste0(substring(col, 1, 7), "4D")
    } else {
      # Named color or other format -> convert via col2rgb
      rgb_vals <- tryCatch(
        grDevices::col2rgb(col),
        error = function(e) return(col)
      )
      if (is.character(rgb_vals)) return(rgb_vals)
      sprintf("#%02X%02X%02X4D", rgb_vals[1], rgb_vals[2], rgb_vals[3])
    }
  }, character(1), USE.NAMES = FALSE)
}

#' Styling and finalizing pedigree graph
#' @import data.table
#' @keywords internal
get_highlight_ids <- function(ped, highlight, trace) {
  focal_ids <- NULL
  relative_ids <- NULL
  trace_edges <- data.table(parent = character(0), child = character(0))
  
  # Default colors
  colors <- list(
    focal_frame = "#9c27b0",
    focal_fill = NULL,
    rel_frame = "#ab47bc",
    rel_fill = NULL
  )
  
  if (!is.null(highlight)) {
    if (is.list(highlight)) {
      focal_ids <- highlight$ids
      if (!is.null(highlight$frame.color)) colors$focal_frame <- highlight$frame.color
      if (!is.null(highlight$color)) colors$focal_fill <- highlight$color
      if (!is.null(highlight$rel.frame.color)) colors$rel_frame <- highlight$rel.frame.color
      if (!is.null(highlight$rel.color)) colors$rel_fill <- highlight$rel.color
    } else {
      focal_ids <- highlight
    }
    
    all_inds <- ped$Ind
    focal_ids <- intersect(focal_ids, all_inds)
    
    if (!isFALSE(trace) && length(focal_ids) > 0) {
      trace_dir <- if (isTRUE(trace)) "all" else trace
      relatives_ped <- suppressWarnings(tidyped(ped, cand = focal_ids, trace = trace_dir))
      relative_ids <- setdiff(unique(relatives_ped$Ind), focal_ids)
      
      # Build trace_edges: parent-child pairs that are ON the trace path.
      # For trace="all", we must call tidyped separately for "up" and "down"
      # to avoid cross-path edges (e.g. N is X's ancestor via Q->T->U,
      # but N->Z1/Z2 is NOT on X's trace path even though both N and Z1/Z2
      # are highlighted).
      if (trace_dir == "all") {
        ped_up <- suppressWarnings(tidyped(ped, cand = focal_ids, trace = "up"))
        ped_down <- suppressWarnings(tidyped(ped, cand = focal_ids, trace = "down"))
        trace_edges <- rbind(
          ped_up[!is.na(Sire), .(parent = Sire, child = Ind)],
          ped_up[!is.na(Dam), .(parent = Dam, child = Ind)],
          ped_down[!is.na(Sire), .(parent = Sire, child = Ind)],
          ped_down[!is.na(Dam), .(parent = Dam, child = Ind)]
        )
      } else {
        trace_edges <- rbind(
          relatives_ped[!is.na(Sire), .(parent = Sire, child = Ind)],
          relatives_ped[!is.na(Dam), .(parent = Dam, child = Ind)]
        )
      }
      trace_edges <- unique(trace_edges)
    }
  }
  
  list(focal = focal_ids, relatives = relative_ids, all_ids = c(focal_ids, relative_ids),
       colors = colors, trace_edges = trace_edges)
}

#' Apply node styles (color, shape, highlighting)
#' @import data.table
#' @keywords internal
apply_node_styles <- function(ped_node, highlight_info) {
  shape = frame.color = color = size = label.color = frame.width = label.font = highlighted = NULL
  
  # Default styles (set all at once)
  ped_node[, `:=`(
    shape = "circle", 
    frame.color = "#7fae59", 
    color = "#9cb383", 
    size = 15, 
    label.color = "#0d0312", 
    frame.width = 0.2, 
    label.font = 1, 
    highlighted = FALSE
  )]
  
  # Using sub-selections is faster than multiple := lines
  ped_node[nodetype == "compact", shape := "square"]
  ped_node[sex == "male", `:=`(frame.color = "#119ecc", color = "#119ecc")]
  ped_node[sex == "female", `:=`(frame.color = "#f4b131", color = "#f4b131")]
  ped_node[is.na(sex) | sex == "unknown", `:=`(frame.color = "#9cb383", color = "#9cb383")]
  
  # Virtual nodes: circle with tiny size and transparency to fix edge gaps
  ped_node[nodetype == "virtual", `:=`(
    shape = "circle", 
    label = "", 
    size = 0.001, 
    frame.width = 0, 
    color = "#FFFFFF00", 
    frame.color = "#FFFFFF00"
  )]
  
  # Apply highlighting
  h_ids <- highlight_info$all_ids
  if (length(h_ids) > 0) {
    focal <- highlight_info$focal
    relatives <- highlight_info$relatives
    colors <- highlight_info$colors
    
    # Use binary search / keys for matching if many ids
    if (length(relatives) > 0) {
      ped_node[label %in% relatives & nodetype == "real", highlighted := TRUE]
    }
    if (length(focal) > 0) {
      ped_node[label %in% focal & nodetype == "real", highlighted := TRUE]
    }
    
    h_familynums <- ped_node[highlighted == TRUE, unique(familynum)]
    ped_node[id %in% h_familynums & nodetype == "virtual", highlighted := TRUE]
    
    # Batch update non-highlighted
    ped_node[highlighted == FALSE & nodetype %in% c("real", "compact"), `:=`(
      color = fade_cols(color), 
      frame.color = fade_cols(frame.color), 
      label.color = fade_cols(label.color)
    )]
    
    # Highlighted styles
    if (length(relatives) > 0) {
      ped_node[label %in% relatives & nodetype == "real", `:=`(
        frame.color = colors$rel_frame, 
        frame.width = 0.5,
        color = if (!is.null(colors$rel_fill)) colors$rel_fill else color
      )]
    }
    
    if (length(focal) > 0) {
      ped_node[label %in% focal & nodetype == "real", `:=`(
        frame.color = colors$focal_frame, 
        frame.width = 1, 
        label.font = 2,
        color = if (!is.null(colors$focal_fill)) colors$focal_fill else color
      )]
    }
  }
  return(ped_node)
}

#' Finalize graph and reindex IDs
#' @import data.table
#' @keywords internal
finalize_graph <- function(ped_node, ped_edge, highlight_info, trace, showf) {
  old_ids <- ped_node$id
  ped_node[, id := seq_len(.N)]
  ped_edge[, from := ped_node$id[match(from, old_ids)]]
  ped_edge[, to := ped_node$id[match(to, old_ids)]]
  
  real_max <- max(ped_node[nodetype %in% c("real", "compact")]$id, na.rm = TRUE)
  
  tonodecolor = i.color = i.highlighted = from_highlighted = NULL
  ped_edge[ped_node, ":="(tonodecolor = i.color, to_highlighted = i.highlighted), on = .(to = id)]
  ped_edge[ped_node, from_highlighted := i.highlighted, on = .(from = id)]
  
  h_ids <- highlight_info$all_ids
  has_trace <- !isFALSE(trace) && length(highlight_info$relatives) > 0
  trace_edges <- highlight_info$trace_edges
  
  # Default: edges from family nodes to parents follow the parent node color
  ped_edge[from > real_max, color := tonodecolor]
  
  # If highlighting is active and family node is not highlighted, fade the edge
  if (length(h_ids) > 0) {
    ped_edge[from > real_max & from_highlighted == FALSE, color := fade_cols(tonodecolor)]
  }
  
  if (length(h_ids) > 0 && has_trace) {
    # Use trace_edges to precisely control which family->parent edges are highlighted.
    # A family->parent edge should only be highlighted if there exists at least one
    # child of that family whose parent-child relationship is on the trace path.
    if (!is.null(trace_edges) && nrow(trace_edges) > 0) {
      # Get the label of the 'to' node (the parent)
      ped_edge[ped_node, to_label := i.label, on = .(to = id)]
      # Get the label of the 'from' node (the family node's label = familylabel)
      ped_edge[ped_node, from_label := i.label, on = .(from = id)]
      
      # For each individual->family edge (from <= real_max), record child-family mapping
      child_to_family <- ped_edge[from <= real_max, .(child_id = from, family_id = to)]
      child_to_family[ped_node, child_label := i.label, on = .(child_id = id)]
      
      # For family->parent edges that are currently highlighted (from_highlighted == TRUE),
      # check if the parent-child relationship is actually on the trace path
      fam_parent_idx <- which(ped_edge$from > real_max & ped_edge$from_highlighted == TRUE)
      
      for (i in fam_parent_idx) {
        fam_id <- ped_edge$from[i]
        parent_label <- ped_edge$to_label[i]
        
        # Find children connected to this family node
        children_labels <- child_to_family[family_id == fam_id, child_label]
        
        # Check if any child-parent pair is in trace_edges
        on_path <- any(trace_edges[parent == parent_label & child %in% children_labels, .N] > 0)
        
        if (!on_path) {
          # This edge is NOT on the trace path - fade it
          ped_edge[i, color := fade_cols(ped_edge$tonodecolor[i])]
        }
      }
      
      # Clean up temporary columns
      ped_edge[, c("to_label", "from_label") := NULL]
    }
    
    # For individual->family edges: highlight only if this individual
    # appears as a child in trace_edges (i.e., the upward connection
    # to its parents is on the trace path).
    ped_edge[ped_node, from_label := i.label, on = .(from = id)]
    ped_edge[
      from <= real_max & from_highlighted == TRUE,
      color := fifelse(from_label %in% trace_edges$child, "#333333", "#3333334D")
    ]
    ped_edge[, from_label := NULL]
  } else if (length(h_ids) == 0) {
    # No highlighting: edges already follow parent node color
  }
  # When highlighting without trace (has_trace == FALSE), keep individual->family edges faded
  
  ped_edge[from <= real_max, ":="(curved = 0, arrow.size = 0, arrow.width = 0, arrow.mode = 0)]
  
  new_names_edge <- c(c("from", "to"), setdiff(colnames(ped_edge), c("from", "to", "tonodecolor")))
  ped_edge <- ped_edge[, ..new_names_edge][order(from, to)]
  
  new_names_node <- c("id", setdiff(colnames(ped_node), "id"))
  ped_node <- ped_node[, ..new_names_node][order(layer, id)]
  
  # Ensure Ind column exists for layout matching (clean label before modification)
  ped_node[, Ind := label]
  
  if (showf && "f" %in% colnames(ped_node)) {
    ped_node[nodetype %in% c("real", "compact") & !is.na(f) & f > 0, 
             label := paste0(label, "\n[", round(f, 4), "]")]
  }
  list(node = ped_node, edge = ped_edge)
}
