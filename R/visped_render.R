#' Left margin (in normalized [0,1] coords) reserved for generation labels
#' @keywords internal
GENLAB_MARGIN <- 0.05

#' Generation label color
#' @keywords internal
GENLAB_COLOR <- "#888888"

#' Render pedigree graph using Two-Pass strategy
#' @importFrom igraph V E plot.igraph vertex_attr edge_attr
#' @importFrom utils modifyList
#' @importFrom stats median
#' @keywords internal
plot_ped_igraph <- function(g, l, node_size, gen_info = NULL, genlab = FALSE, ...) {
  # ============================================================================
  # SCALING LOGIC FOR RENDERING
  # ============================================================================
  scaling_ref <- max(node_size, 1.5)
  
  # Access attributes as lists/vectors ONCE to avoid repeated V(g)/E(g) lookups
  # Using vertex_attr/edge_attr is generally more stable for large graphs
  v_attrs <- igraph::vertex_attr(g)
  e_attrs <- igraph::edge_attr(g)
  
  # Scale frame width
  v_frame_width_base <- if (!is.null(v_attrs$frame.width)) v_attrs$frame.width else 0.2
  lwd_scale <- if (scaling_ref < 2.5) scaling_ref / 2.5 else 1
  v_frame_width <- v_frame_width_base * lwd_scale
  
  # Scale edge width
  e_width_base <- if (!is.null(e_attrs$width)) e_attrs$width else 1
  e_width <- pmax(e_width_base * (scaling_ref * 0.15), 0.001)
  
  # Scale arrows
  e_arrow_size_base <- if (!is.null(e_attrs$arrow.size)) e_attrs$arrow.size else 1
  calc_arrow_size <- e_arrow_size_base * (scaling_ref * 0.005)
  
  e_arrow_size <- calc_arrow_size
  e_arrow_width <- if (!is.null(e_attrs$arrow.width)) e_attrs$arrow.width else 1
  e_arrow_mode <- if (!is.null(e_attrs$arrow.mode)) e_attrs$arrow.mode else 2
  
  if (median(calc_arrow_size, na.rm = TRUE) < 0.005) {
    e_arrow_size <- 0
    e_arrow_width <- 0
    e_arrow_mode <- 0
  } else {
    e_arrow_width <- e_arrow_width * (scaling_ref * 0.015)
  }
  
  user_args <- list(...)
  margin <- max(node_size / 100, 0.02)
  genlab_space <- if (genlab && !is.null(gen_info) && nrow(gen_info) > 0) GENLAB_MARGIN else 0
  
  # PASS 1: Draw EDGES ONLY
  # Instead of modifying g, we pass visual attributes directly to plot.igraph
  plot_args_edges <- list(
    x = g, rescale = FALSE,
    xlim = c(-margin - genlab_space, 1 + margin), ylim = c(1 + margin, -margin),
    layout = l, asp = 0, add = FALSE,
    # Edge styles
    edge.width = e_width,
    edge.color = if (!is.null(e_attrs$color)) e_attrs$color else "#333333",
    edge.curved = if (!is.null(e_attrs$curved)) e_attrs$curved else 0.1,
    edge.arrow.size = e_arrow_size,
    edge.arrow.width = e_arrow_width,
    edge.arrow.mode = e_arrow_mode,
    # Hide nodes
    vertex.size = 0,
    vertex.label = NA,
    vertex.color = NA,
    vertex.frame.color = NA,
    vertex.frame.width = 0
  )
  
  if (length(user_args) > 0) {
    # Filter out user vertex args from edge pass
    safe_user_args <- user_args[!grepl("^(vertex|label)\\.", names(user_args))]
    plot_args_edges <- utils::modifyList(plot_args_edges, safe_user_args)
    # Re-enforce node hiding for pass 1
    plot_args_edges$vertex.size <- 0
    plot_args_edges$vertex.label <- NA
  }
  
  suppressWarnings(do.call(igraph::plot.igraph, plot_args_edges))
  
  # PASS 2: Draw NODES AND LABELS ONLY
  # We reuse the same graph object but hide edges
  plot_args_nodes <- list(
    x = g, rescale = FALSE,
    xlim = c(-margin - genlab_space, 1 + margin), ylim = c(1 + margin, -margin),
    layout = l, asp = 0, add = TRUE,
    # Node styles
    vertex.size = if (!is.null(v_attrs$size)) v_attrs$size else 15,
    vertex.shape = if (!is.null(v_attrs$shape)) v_attrs$shape else "circle",
    vertex.color = if (!is.null(v_attrs$color)) v_attrs$color else "#9cb383",
    vertex.frame.color = if (!is.null(v_attrs$frame.color)) v_attrs$frame.color else "#7fae59",
    vertex.frame.width = v_frame_width,
    vertex.label = if (!is.null(v_attrs$label)) v_attrs$label else v_attrs$name,
    vertex.label.cex = if (!is.null(v_attrs$label.cex)) v_attrs$label.cex else 1,
    vertex.label.color = if (!is.null(v_attrs$label.color)) v_attrs$label.color else "#0d0312",
    vertex.label.font = if (!is.null(v_attrs$label.font)) v_attrs$label.font else 1,
    # Hide edges
    edge.color = "#FFFFFF00",
    edge.width = 0,
    edge.arrow.size = 0,
    edge.arrow.mode = 0
  )
  
  if (length(user_args) > 0) {
    plot_args_nodes <- utils::modifyList(plot_args_nodes, user_args)
  }
  
  suppressWarnings(do.call(igraph::plot.igraph, plot_args_nodes))

  # PASS 3: Draw generation labels on the left margin
  if (genlab && !is.null(gen_info) && nrow(gen_info) > 0) {
    gen_label_cex <- max(0.5, min(1.5, scaling_ref * 0.3))
    gen_label_x <- -margin - genlab_space * 0.5
    for (i in seq_len(nrow(gen_info))) {
      graphics::text(
        x = gen_label_x,
        y = gen_info$y[i],
        labels = paste0("G", gen_info$gen[i]),
        cex = gen_label_cex,
        col = GENLAB_COLOR,
        font = 2,
        adj = c(0.5, 0.5)
      )
    }
  }
}

