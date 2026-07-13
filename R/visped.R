#' Normalize custom labels before pedigree tidying
#' @keywords internal
prepare_visped_labels <- function(ped, labelvar) {
  if (is.null(labelvar)) {
    return(list(ped = ped, label_column = NULL))
  }

  if (!is.character(labelvar)) {
    stop(
      "'labelvar' must be NULL, a single non-empty column name, ",
      "or a character vector with one label per pedigree row."
    )
  }

  if (length(labelvar) == 1L) {
    if (is.na(labelvar) || labelvar == "") {
      stop(
        "'labelvar' must be NULL, a single non-empty column name, ",
        "or a character vector with one label per pedigree row."
      )
    }
    return(list(ped = ped, label_column = labelvar))
  }

  if (!is.data.frame(ped)) {
    stop("A 'labelvar' vector requires 'ped' to be a data.frame or data.table.")
  }

  if (length(labelvar) != nrow(ped)) {
    stop(
      sprintf(
        "A 'labelvar' vector must have one label per pedigree row (%d), not %d.",
        nrow(ped),
        length(labelvar)
      )
    )
  }

  ped <- copy(ped)
  label_column <- ".visped_label"
  while (label_column %in% names(ped)) {
    label_column <- paste0(label_column, "_")
  }
  set(ped, j = label_column, value = labelvar)

  list(ped = ped, label_column = label_column)
}

#' Normalize generation-label display settings
#' @keywords internal
prepare_visped_genlab <- function(genlab) {
  if (is.logical(genlab) && length(genlab) == 1L && !is.na(genlab)) {
    return(list(show = genlab, labels = NULL))
  }

  if (is.character(genlab) && length(genlab) > 0L) {
    if (anyNA(genlab) || any(!nzchar(genlab))) {
      stop("Custom 'genlab' labels must not be NA or empty strings.")
    }

    label_names <- names(genlab)
    if (!is.null(label_names) && any(!is.na(label_names) & nzchar(label_names))) {
      stop("Custom 'genlab' labels must be an unnamed character vector.")
    }

    return(list(show = TRUE, labels = unname(genlab)))
  }

  stop(
    "'genlab' must be TRUE, FALSE, or a non-empty unnamed character vector."
  )
}

#' Attach display labels to the final plotted generation layers
#' @keywords internal
prepare_visped_geninfo <- function(gen_info, labels = NULL) {
  gen_info <- copy(gen_info)

  if (!is.null(labels) && length(labels) != nrow(gen_info)) {
    stop(sprintf(
      "Custom 'genlab' labels must have one label for each displayed generation (%d), not %d.",
      nrow(gen_info),
      length(labels)
    ))
  }

  gen_info[, label := if (is.null(labels)) paste0("G", gen) else labels]
  return(gen_info[])
}

#' Visualize a tidy pedigree
#'
#' \code{visped} function draws a graph of a full or compact pedigree.
#'
#' This function takes a pedigree tidied by the \code{\link{tidyped}} function and outputs a hierarchical graph for all individuals in the pedigree. The graph can be shown on the default graphic device or saved as a PDF or SVG file. The vector output is legible and avoids overlapping labels. It is especially useful when the number of individuals is large and individual labels are long.
#'
#' Rendering is performed using a Two-Pass strategy: edges are drawn first to ensure center-to-center connectivity, followed by nodes and labels. This ensures perfect visual alignment in high-resolution vector outputs. The function also supports real-time ancestry and descendant highlighting.
#'
#' This function can draw the graph of a very large pedigree (> 10,000 individuals per generation) by compacting full-sib individuals. It is highly effective for aquatic animal pedigrees, which usually include many full-sib families per generation in nucleus breeding populations. The outline of a pedigree without individual labels is still shown if the width of a pedigree graph exceeds the maximum width (500 inches) of the PDF file.
#'
#' By default, \code{shapeby = "sex"} encodes individual sex with node shape:
#' females are circles, males are squares, individuals of unknown sex are
#' diamonds, and monoecious individuals are hexagons. Set
#' \code{shapeby = "role"} to use the legacy role-based scheme, where circles
#' represent individual records and green-grey rectangles represent compact
#' full-sib family summaries. Compact family summaries remain green-grey
#' rectangles in both modes. In \code{shapeby = "role"} mode, real individuals
#' remain circles, but their fill colors still encode sex. Dark sky blue
#' indicates males, dark goldenrod indicates females, teal indicates monoecious
#' individuals, and neutral grey indicates unknown sex. Purple borders indicate
#' highlighted individuals rather than sex.
#'
#' @param ped A \code{tidyped} object (which inherits from \code{data.table}). It is recommended that the pedigree is tidied and pruned by candidates using the \code{\link{tidyped}} function with the non-null parameter \code{cand}.
#' @param compact A logical value indicating whether terminal, non-parent
#'   full-sib individuals will be replaced by one family display node. For
#'   example, 100 eligible full siblings are shown as \code{"FS×100"} when
#'   \code{compact = TRUE}. Highlighted individuals are not compacted, and the
#'   underlying pedigree is unchanged. The default value is FALSE.
#' @param outline A logical value indicating whether shapes without labels will be shown. A graph of the pedigree without individual labels is shown when setting \code{outline = TRUE}. This is useful for viewing the pedigree outline and identifying immigrant individuals in each generation when the graph width exceeds the maximum PDF width (500 inches). The default value is FALSE.
#' @param cex NULL or a numeric value changing the size of individual labels shown in the graph. \emph{cex} is an abbreviation for 'character expansion factor'. The \code{visped} function will attempt to estimate (\code{cex=NULL}) the appropriate cex value and report it in the messages. Based on the reported cex from a previous run, this parameter should be increased if labels are wider than their shapes in the PDF; conversely, it should be decreased if labels are narrower than their shapes. The default value is NULL.
#' @param showgraph A logical value indicating whether a plot will be shown in the default graphic device (e.g., the Plots panel in RStudio). This is useful for quick viewing without opening a PDF file. However, the graph on the default device may not be legible (e.g., overlapping labels or aliasing lines) due to size restrictions. It is recommended to set \code{showgraph = FALSE} for large pedigrees. The default value is TRUE.
#' @param file NULL or a character value specifying whether the pedigree graph will be saved as a vector file. Files ending in '.svg' are written with the SVG device; all other file names use PDF output. The vector output is legible and avoids overlapping labels even with many individuals or long labels. It is recommended to save the pedigree graph as a vector file. The default value is NULL.
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
#' @param labelvar NULL, a single non-empty character value naming a column in
#'   \code{ped}, or a character vector with one displayed label per input
#'   pedigree row. Row-aligned vectors are attached before any internal
#'   tidying or reordering. Compact full-sib family nodes still show family
#'   size when \code{compact = TRUE}. Missing or empty values fall back to
#'   individual IDs. The default value is NULL.
#' @param shapeby Character string controlling what node shapes encode.
#'   \code{"sex"} (default) uses circles for females, squares for males,
#'   diamonds for unknown sex, and hexagons for monoecious individuals.
#'   \code{"role"} distinguishes individual records from compact family
#'   summaries by shape while keeping sex-based fill colors for real
#'   individuals. Compact family summaries are rectangles in both modes.
#' @param pagewidth A numeric value specifying the width of the PDF file in inches. This controls the horizontal scaling of the layout. The default value is 200.
#' @param symbolsize A numeric value specifying the scaling factor for node size relative to the label size. Values greater than 1 increase the node size (adding padding around the label), while values less than 1 decrease it. This is useful for fine-tuning the whitespace and legibility of dense graphs. The default value is 1.
#' @param maxiter An integer specifying the maximum number of iterations for the Sugiyama layout algorithm to minimize edge crossings. Higher values (e.g., 2000 or 5000) may result in fewer crossed lines for complex pedigrees but will increase computation time. The default value is 1000.
#' @param genlab FALSE to omit generation labels; TRUE to draw the default
#'   labels G1, G2, ...; or an unnamed, non-empty character vector of custom
#'   labels. Custom labels are assigned from top to bottom to the final
#'   displayed generation layers and must contain exactly one label per layer.
#'   The default value is FALSE.
#' @param genlabcex NULL or a numeric value controlling the size of displayed
#'   generation labels. If \code{NULL}, \code{visped()} uses an automatic size
#'   based on node scaling. Set a larger value to keep generation labels
#'   readable in deep pedigrees. The default value is NULL.
#' @param ... Additional arguments passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @return The function mainly produces a plot on the current graphics device and/or a vector file. It invisibly returns a list containing the graph object, layout coordinates, and node sizes.
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
#' # Use application-specific labels for the displayed generations.
#' # Custom labels are assigned from top to bottom.
#' generation_labels <- paste("Generation", sort(unique(simple_ped_tidy$Gen)))
#' visped(simple_ped_tidy,
#'        genlab = generation_labels,
#'        cex = 0.25,
#'        symbolsize = 5.5)
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
#' # Use the role-based legacy symbol scheme
#' visped(simple_ped_tidy,
#'        shapeby = "role",
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
#' # Saving pedigree as SVG with larger labels and tighter symbols
#' visped(simple_ped_tidy,
#'        cex=0.8,
#'        symbolsize=1,
#'        pagewidth = 100,
#'        file = tempfile(fileext = ".svg"))
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
#' @importFrom grDevices pdf svg dev.off dev.cur
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
  labelvar = NULL,
  shapeby = c("sex", "role"),
  pagewidth = 200,
  symbolsize = 1,
  maxiter = 1000,
  genlab = FALSE,
  genlabcex = NULL,
  ...
) {
  label_info <- prepare_visped_labels(ped, labelvar)
  ped <- label_info$ped
  label_column <- label_info$label_column
  genlab_info <- prepare_visped_genlab(genlab)

  # Automatically convert raw data to tidyped object if needed.
  # If the object already looks like a tidyped pedigree but only lost its class,
  # restore/validate it instead of rebuilding from raw input.
  tidyped_core <- c("Ind", "Sire", "Dam", "Sex", "Gen", "IndNum", "SireNum", "DamNum")
  looks_tidyped <- is.data.frame(ped) && all(tidyped_core %in% colnames(ped))

  if (!is_tidyped(ped) && !looks_tidyped) {
    ped <- tidyped(ped, addgen = TRUE, addnum = TRUE)
  }

  ped <- validate_tidyped(ped)
  shapeby <- match.arg(shapeby)

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

  if (!is.null(label_column) && !label_column %in% names(ped)) {
    stop(sprintf("Column '%s' specified by 'labelvar' was not found in the pedigree.", label_column))
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
          labelvar = label_column,
          shapeby = shapeby,
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
      labelvar = label_column,
      shapeby = shapeby,
      pagewidth = pagewidth,
      symbolsize = symbolsize,
      maxiter = maxiter,
      ...
    )
  }

  g <- graph_data$g
  l <- graph_data$layout
  node_size <- graph_data$node_size
  if (genlab_info$show) {
    graph_data$gen_info <- prepare_visped_geninfo(
      graph_data$gen_info,
      genlab_info$labels
    )
  }

  #===Draw the pedigree================================================================
  if (showgraph) {
    plot_ped_igraph(
      g,
      l,
      node_size,
      gen_info = graph_data$gen_info,
      genlab = genlab_info$show,
      genlabcex = genlabcex,
      ...
    )
  }

  if (!is.null(file)) {
    file_ext <- tolower(tools::file_ext(file))
    output_device <- if (identical(file_ext, "svg")) grDevices::svg else grDevices::pdf
    output_device(
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
      genlab = genlab_info$show,
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
    current_genlabcex <- if (
      genlab_info$show &&
        !is.null(graph_data$gen_info) &&
        nrow(graph_data$gen_info) > 0L
    ) {
      resolve_visped_genlabcex(node_size, genlabcex)
    } else {
      NULL
    }
    message(paste(
      "Label cex: ",
      current_cex,
      ". Symbol size: ",
      symbolsize,
      if (!is.null(current_genlabcex)) {
        paste0(". Generation label cex: ", current_genlabcex)
      },
      if (!is.null(current_genlabcex)) {
        ". Adjust 'cex', 'symbolsize', and 'genlabcex' if labels are too large or small."
      } else {
        ". Adjust 'cex' and 'symbolsize' if labels are too large or small."
      },
      sep = ""
    ))
  }

  if ((showgraph || !is.null(file)) && is.null(file)) {
    message("Tip: Use 'file' to save as a legible vector PDF or SVG.")
  }

  if (showf && "f" %in% colnames(ped) && any(ped$f == 0, na.rm = TRUE)) {
    message("Note: Inbreeding coefficients of 0 are not shown in the graph.")
  }

  graph_data$gen_info <- NULL
  invisible(graph_data)
}
