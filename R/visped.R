#' Visualize a tidy pedigree
#'
#' \code{visped} function draws a graph of a full or compact pedigree.
#'
#' This function takes a pedigree tidied by the \code{\link{tidyped}} function and outputs a hierarchical graph for all individuals in the pedigree. The graph can be shown on the default graphic device or saved as a PDF file. The PDF output is a legible vector drawing that is legible and avoids overlapping labels. It is especially useful when the number of individuals is large and individual labels are long.
#'
#' Rendering is performed using a Two-Pass strategy: edges are drawn first to ensure center-to-center connectivity, followed by nodes and labels. This ensures perfect visual alignment in high-resolution vector outputs. The function also supports real-time ancestry and descendant highlighting.
#'
#' This function can draw the graph of a very large pedigree (> 10,000 individuals per generation) by compacting full-sib individuals. It is highly effective for aquatic animal pedigrees, which usually include many full-sib families per generation in nucleus breeding populations. The outline of a pedigree without individual labels is still shown if the width of a pedigree graph exceeds the maximum width (500 inches) of the PDF file.
#'
#' In the graph, two shapes and four colors are used. Circles represent individuals, and squares represent families. Dark sky blue indicates males, dark goldenrod indicates females, purple indicates monoecious individuals (common in plant breeding, where the same individual serves as both male and female parent), and dark olive green indicates unknown sex. For example, a dark sky blue circle represents a male individual; a dark goldenrod square represents all female individuals in a full-sib family when \code{compact = TRUE}.
#'
#' @param ped A \code{tidyped} object (which inherits from \code{data.table}). It is recommended that the pedigree is tidied and pruned by candidates using the \code{\link{tidyped}} function with the non-null parameter \code{cand}.
#' @param compact A logical value indicating whether IDs of full-sib individuals in one generation will be removed and replaced with the number of full-sib individuals. For example, if there are 100 full-sib individuals in one generation, they will be replaced with a single label "100" when \code{compact = TRUE}. The default value is FALSE.
#' @param outline A logical value indicating whether shapes without labels will be shown. A graph of the pedigree without individual labels is shown when setting \code{outline = TRUE}. This is useful for viewing the pedigree outline and identifying immigrant individuals in each generation when the graph width exceeds the maximum PDF width (500 inches). The default value is FALSE.
#' @param cex NULL or a numeric value changing the size of individual labels shown in the graph. \emph{cex} is an abbreviation for 'character expansion factor'. The \code{visped} function will attempt to estimate (\code{cex=NULL}) the appropriate cex value and report it in the messages. Based on the reported cex from a previous run, this parameter should be increased if labels are wider than their shapes in the PDF; conversely, it should be decreased if labels are narrower than their shapes. The default value is NULL.
#' @param showgraph A logical value indicating whether a plot will be shown in the default graphic device (e.g., the Plots panel in RStudio). This is useful for quick viewing without opening a PDF file. However, the graph on the default device may not be legible (e.g., overlapping labels or aliasing lines) due to size restrictions. It is recommended to set \code{showgraph = FALSE} for large pedigrees. The default value is TRUE.
#' @param file NULL or a character value specifying whether the pedigree graph will be saved as a PDF file. The PDF output is a legible vector drawing where labels do not overlap, even with many individuals or long labels. It is recommended to save the pedigree graph as a PDF file. The default value is NULL.
#' @param highlight NULL, a character vector of individual IDs, or a list specifying individuals to highlight. If a character vector is provided, individuals will be highlighted with a purple border while preserving their sex-based fill color. If a list is provided, it should contain:
#' \itemize{
#'   \item \code{ids}: (required) character vector of individual IDs to highlight.
#'   \item \code{frame.color}: (optional) hex color for the border of focal individuals.
#'   \item \code{color}: (optional) hex color for the fill of focal individuals.
#'   \item \code{rel.frame.color}: (optional) hex color for the border of relatives (used when \code{trace} is not NULL).
#'   \item \code{rel.color}: (optional) hex color for the fill of relatives (used when \code{trace} is not NULL).
#' }
#' For example: \code{c("A", "B")} or \code{list(ids = c("A", "B"), frame.color = "#9c27b0")}. The function will check if the specified individuals exist in the pedigree and issue a warning for any missing IDs. The default value is NULL.
#' @param trace A logical value or a character string. If TRUE, all ancestors and descendants of the individuals specified in \code{highlight} will be highlighted. If a character string, it specifies the tracing direction: "\strong{up}" (ancestors), "\strong{down}" (descendants), or "\strong{all}" (union of ancestors and descendants). This is useful for focusing on specific families within a large pedigree. The default value is FALSE.
#' @param showf A logical value indicating whether inbreeding coefficients will be shown in the graph. If \code{showf = TRUE} and the column \strong{f} is missing, \code{visped()} will try to compute it automatically with \code{\link{inbreed}} on a structurally complete pedigree. If automatic computation is not possible, a warning is issued and labels are drawn without \strong{f}. The default value is FALSE.
#' @param pagewidth A numeric value specifying the width of the PDF file in inches. This controls the horizontal scaling of the layout. The default value is 200.
#' @param symbolsize A numeric value specifying the scaling factor for node size relative to the label size. Values greater than 1 increase the node size (adding padding around the label), while values less than 1 decrease it. This is useful for fine-tuning the whitespace and legibility of dense graphs. The default value is 1.
#' @param maxiter An integer specifying the maximum number of iterations for the Sugiyama layout algorithm to minimize edge crossings. Higher values (e.g., 2000 or 5000) may result in fewer crossed lines for complex pedigrees but will increase computation time. The default value is 1000.
#' @param genlab A logical value indicating whether generation labels (G1, G2, ...) will be drawn on the left margin of the pedigree graph. This helps identify the generation of each row of nodes, especially in deep pedigrees with many generations. The default value is FALSE.
#' @param genlabcex NULL or a numeric value controlling the size of generation labels shown when \code{genlab = TRUE}. If \code{NULL}, \code{visped()} uses an automatic size based on node scaling. Set a larger value to keep generation labels readable in deep pedigrees. The default value is NULL.
#' @param ... Additional arguments passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @return The function mainly produces a plot on the current graphics device and/or a PDF file. It invisibly returns a list containing the graph object, layout coordinates, and node sizes.
#'
#' @seealso
#' \code{\link{tidyped}} for tidying pedigree data (required input)
#' \code{\link{vismat}} for visualizing relationship matrices as heatmaps
#' \code{\link{pedmat}} for computing relationship matrices
#' \code{\link{splitped}} for splitting pedigree into connected components
#' \code{\link[igraph:plot.igraph]{plot.igraph}} underlying plotting function
#'
#' @note Isolated individuals (those with no parents and no progeny, assigned Gen 0) are automatically filtered out and not shown in the plot. A message will be issued if any such individuals are removed.
#'
#' @examples
#' library(visPedigree)
#' library(data.table)
#' # Drawing a simple pedigree
#' simple_ped_tidy <- tidyped(simple_ped)
#' visped(simple_ped_tidy, 
#'        cex=0.25, 
#'        symbolsize=5.5)
#'
#' # Highlighting an individual and its ancestors and descendants
#' visped(simple_ped_tidy, 
#'        highlight = "J5X804", 
#'        trace = "all", 
#'        cex=0.25, 
#'        symbolsize=5.5)
#'
#' # Showing inbreeding coefficients in the graph
#' simple_ped_tidy_inbreed <- tidyped(simple_ped, inbreed = TRUE)
#' visped(simple_ped_tidy_inbreed,
#'        showf = TRUE, 
#'        cex=0.25, 
#'        symbolsize=5.5)
#'
#' # visped() will automatically compute inbreeding coefficients if 'f' is missing
#' visped(simple_ped_tidy,
#'        showf = TRUE,
#'        cex=0.25,
#'        symbolsize=5.5)
#'
#' # Adjusting page width and symbol size for better layout
#' # Increase pagewidth to spread nodes horizontally in the pdf file
#' # Increase symbolsize for more padding around individual labels
#' visped(simple_ped_tidy, 
#'        cex=0.25, 
#'        symbolsize=5.5, 
#'        pagewidth = 100, 
#'        file = tempfile(fileext = ".pdf"))
#'
#' # Highlighting multiple individuals with custom colors
#' visped(simple_ped_tidy,
#'        highlight = list(ids = c("J3Y620", "J1X971"),
#'                         frame.color = "#4caf50",
#'                         color = "#81c784"),
#'        cex=0.25,
#'        symbolsize=5.5)
#'
#' # Handling large pedigrees: Saving to PDF is recommended for legibility
#' # The 'trace' and 'tracegen' parameters in tidyped() help prune the graph
#' cand_labels <- big_family_size_ped[(Year == 2007) & (substr(Ind,1,2) == "G8"), Ind]
#' \donttest{
#' big_ped_tidy <- tidyped(big_family_size_ped, 
#'                         cand = cand_labels, 
#'                         trace = "up", 
#'                         tracegen = 2)
#' # Use compact = TRUE for large families
#' visped(big_ped_tidy, 
#'        compact = TRUE, 
#'        cex=0.08, 
#'        symbolsize=5.5, 
#'        file = tempfile(fileext = ".pdf"))
#'
#' # Use outline = TRUE if individual labels are not required
#' visped(big_ped_tidy, 
#'        compact = TRUE, 
#'        outline = TRUE, 
#'        file = tempfile(fileext = ".pdf"))
#' }
#'
#' @import data.table
#' @import igraph
#' @importFrom grDevices pdf dev.off dev.cur
#' @importFrom graphics strwidth
#' @export
visped <- function(
  ped,
  compact = FALSE,
  outline = FALSE,
  cex = NULL,
  showgraph = TRUE,
  file = NULL,
  highlight = NULL,
  trace = FALSE,
  showf = FALSE,
  pagewidth = 200,
  symbolsize = 1,
  maxiter = 1000,
  genlab = FALSE,
  genlabcex = NULL,
  ...
) {
  # Automatically convert raw data to tidyped object if needed.
  # If the object already looks like a tidyped pedigree but only lost its class,
  # restore/validate it instead of rebuilding from raw input.
  tidyped_core <- c("Ind", "Sire", "Dam", "Sex", "Gen", "IndNum", "SireNum", "DamNum")
  looks_tidyped <- is.data.frame(ped) && all(tidyped_core %in% colnames(ped))

  if (!is_tidyped(ped) && !looks_tidyped) {
    ped <- tidyped(ped, addgen = TRUE, addnum = TRUE)
  }

  ped <- validate_tidyped(ped)

  if (!isTRUE(compact) && !isFALSE(compact)) {
    stop("'compact' must be TRUE or FALSE.")
  }

  if (!isTRUE(outline) && !isFALSE(outline)) {
    stop("'outline' must be TRUE or FALSE.")
  }

  if (!isTRUE(showgraph) && !isFALSE(showgraph)) {
    stop("'showgraph' must be TRUE or FALSE.")
  }

  if (!isTRUE(showf) && !isFALSE(showf)) {
    stop("'showf' must be TRUE or FALSE.")
  }

  if (
    !is.null(cex) &&
      (!is.numeric(cex) || length(cex) != 1 || is.na(cex) || cex <= 0)
  ) {
    stop("'cex' must be NULL or a single positive number.")
  }

  if (
    !is.null(file) && (!is.character(file) || length(file) != 1 || is.na(file))
  ) {
    stop("'file' must be NULL or a single character string.")
  }

  # 1. Ensure at least one output is selected (validated after flags are checked)
  if (!showgraph && is.null(file)) {
    stop(
      "Both 'showgraph' and 'file' are disabled. No output will be generated."
    )
  }

  if (!is.null(highlight) && !is.character(highlight) && !is.list(highlight)) {
    stop("'highlight' must be NULL, a character vector, or a list.")
  }

  if (!isFALSE(trace)) {
    if (isTRUE(trace)) {
      trace <- "all"
    } else if (
      !is.character(trace) ||
        length(trace) != 1 ||
        is.na(trace) ||
        !trace %in% c("up", "down", "all")
    ) {
      stop("'trace' must be TRUE, FALSE, 'up', 'down', or 'all'.")
    }
  }

  if (
    !is.numeric(pagewidth) ||
      length(pagewidth) != 1 ||
      is.na(pagewidth) ||
      pagewidth <= 0
  ) {
    stop("'pagewidth' must be a single positive number.")
  }

  if (
    !is.numeric(symbolsize) ||
      length(symbolsize) != 1 ||
      is.na(symbolsize) ||
      symbolsize <= 0
  ) {
    stop("'symbolsize' must be a single positive number.")
  }

  if (
    !is.numeric(maxiter) ||
      length(maxiter) != 1 ||
      is.na(maxiter) ||
      maxiter <= 0
  ) {
    stop("'maxiter' must be a single positive integer.")
  }
  maxiter <- as.integer(maxiter)

  if (!isTRUE(genlab) && !isFALSE(genlab)) {
    stop("'genlab' must be TRUE or FALSE.")
  }

  if (
    !is.null(genlabcex) &&
      (!is.numeric(genlabcex) || length(genlabcex) != 1 || is.na(genlabcex) || genlabcex <= 0)
  ) {
    stop("'genlabcex' must be NULL or a single positive number.")
  }

  # 2. Sanitize highlight inputs
  if (!is.null(highlight)) {
    if (is.character(highlight)) {
      highlight <- highlight[!is.na(highlight) & highlight != ""]
    } else if (is.list(highlight) && !is.null(highlight[["ids"]])) {
      highlight[["ids"]] <- highlight[["ids"]][
        !is.na(highlight[["ids"]]) & highlight[["ids"]] != ""
      ]
    }
    # If sanitization left it empty, set to NULL
    if (is.character(highlight) && length(highlight) == 0) {
      highlight <- NULL
    }
    if (is.list(highlight) && length(highlight[["ids"]]) == 0) highlight <- NULL
  }

  if (showf && !has_inbreeding(ped)) {
    if (is_complete_pedigree(ped)) {
      ped <- inbreed(ped)
      message(
        "Note: 'showf = TRUE' requested but 'f' column was missing. ",
        "Calculated inbreeding coefficients automatically."
      )
    } else {
      warning(
        "Inbreeding coefficients ('f' column) not found and cannot be ",
        "computed automatically because the pedigree is structurally incomplete. ",
        "Run tidyped(..., inbreed = TRUE) on a complete pedigree first, ",
        "or extract a valid sub-pedigree with tidyped(tp, cand = ids, trace = \"up\")."
      )
      showf <- FALSE
    }
  }

  # Prepare graph data
  # Use a temporary PDF device to calculate text metrics if file output is requested
  # This ensures strwidth/strheight match the target PDF device (e.g. font metrics)
  if (!is.null(file)) {
    tmp_pdf <- tempfile()
    # Pass font family if provided in ... to ensure accurate cex estimation
    dots <- list(...)
    pdf_args <- list(file = tmp_pdf)
    if ("family" %in% names(dots)) {
      pdf_args$family <- dots$family
    }
    do.call(pdf, pdf_args)

    on.exit(
      {
        if (file.exists(tmp_pdf)) unlink(tmp_pdf)
      },
      add = TRUE
    )

    # Calculate graph data using temp device for font metrics
    graph_data <- tryCatch(
      {
        prepare_ped_graph(
          ped = ped,
          compact = compact,
          outline = outline,
          cex = cex,
          highlight = highlight,
          trace = trace,
          showf = showf,
          pagewidth = pagewidth,
          symbolsize = symbolsize,
          maxiter = maxiter,
          ...
        )
      },
      finally = {
        if (dev.cur() > 1) dev.off()
      }
    )
  } else {
    graph_data <- prepare_ped_graph(
      ped = ped,
      compact = compact,
      outline = outline,
      cex = cex,
      highlight = highlight,
      trace = trace,
      showf = showf,
      pagewidth = pagewidth,
      symbolsize = symbolsize,
      maxiter = maxiter,
      ...
    )
  }

  g <- graph_data$g
  l <- graph_data$layout
  node_size <- graph_data$node_size

  #===Draw the pedigree================================================================
  if (showgraph) {
    plot_ped_igraph(
      g,
      l,
      node_size,
      gen_info = graph_data$gen_info,
      genlab = genlab,
      genlabcex = genlabcex,
      ...
    )
  }

  if (!is.null(file)) {
    pdf(
      file = file,
      width = graph_data$canvas_width,
      height = graph_data$canvas_height
    )
    # Ensure the device is closed even if plotting fails.
    on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
    plot_ped_igraph(
      g,
      l,
      node_size,
      gen_info = graph_data$gen_info,
      genlab = genlab,
      genlabcex = genlabcex,
      ...
    )

    # Correct path normalization for the message
    saved_path <- tryCatch(
      normalizePath(file, mustWork = FALSE),
      error = function(e) file
    )
    message(paste("Pedigree saved to: ", saved_path, sep = ""))
  }

  if ((showgraph || !is.null(file)) && !outline) {
    current_cex <- if (is.null(cex)) graph_data$best_cex else cex
    message(paste(
      "Label cex: ",
      current_cex,
      ". Symbol size: ",
      symbolsize,
      ". Adjust 'cex' and 'symbolsize' if labels are too large or small.",
      sep = ""
    ))
  }

  if ((showgraph || !is.null(file)) && is.null(file)) {
    message("Tip: Use 'file' to save as a legible vector PDF.")
  }

  if (showf && "f" %in% colnames(ped) && any(ped$f == 0, na.rm = TRUE)) {
    message("Note: Inbreeding coefficients of 0 are not shown in the graph.")
  }

  graph_data$gen_info <- NULL
  invisible(graph_data)
}
