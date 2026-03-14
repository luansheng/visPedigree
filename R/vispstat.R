#' Visualize Pedigree Statistics
#'
#' \code{vispstat} visualizes statistics from a \code{pedstats} object,
#' including generation intervals and ancestral depth distributions.
#'
#' @param x A \code{pedstats} object returned by \code{\link{pedstats}}.
#' @param type Character. The type of plot to generate:
#' \itemize{
#'   \item \code{"genint"}: Bar chart of generation intervals (Mean ± SD).
#'   \item \code{"ecg"}: Histogram of ancestral depth metrics (ECG, FullGen, or MaxGen).
#' }
#' @param metric Character. Specific metric to plot when \code{type = "ecg"}.
#' Options: \code{"ECG"} (default), \code{"FullGen"}, \code{"MaxGen"}.
#' @param ... Additional arguments passed to \code{\link[lattice]{barchart}} or 
#' \code{\link[lattice]{histogram}}.
#'
#' @return A lattice plot object.
#'
#' @seealso \code{\link{pedstats}}
#' 
#' @examples
#' \dontrun{
#' library(visPedigree)
#' data(simple_ped)
#' 
#' # Add a birth year column for generation interval calculation
#' simple_ped$Year <- sample(2010:2020, nrow(simple_ped), replace = TRUE)
#' tped <- tidyped(simple_ped)
#' 
#' # Calculate statistics
#' stats <- pedstats(tped, timevar = "Year")
#' 
#' # Visualize generation intervals
#' vispstat(stats, type = "genint")
#' 
#' # Visualize ancestral depth (ECG)
#' vispstat(stats, type = "ecg", metric = "ECG")
#' 
#' # Visualize fully traced generations
#' vispstat(stats, type = "ecg", metric = "FullGen")
#' 
#' # Use the plot method
#' plot(stats, type = "ecg", metric = "MaxGen")
#' }
#' 
#' @export
vispstat <- function(x, type = c("genint", "ecg"), metric = "ECG", ...) {
  if (!inherits(x, "pedstats")) {
    stop("x must be a pedstats object")
  }
  
  type <- match.arg(type)
  
  if (type == "genint") {
    if (is.null(x$gen_intervals)) {
      stop("No generation interval data found in pedstats object. ",
           "Ensure the pedigree has a time variable (e.g., 'Year' or 'BirthYear').")
    }
    
    dt <- data.table::copy(x$gen_intervals)
    # Keep only the four classic gametic pathways
    dt <- dt[Pathway %in% c("SS", "SD", "DS", "DD")]
    
    if (nrow(dt) == 0) {
      stop("No pathway-specific generation intervals found.")
    }
    
    # Order pathways consistently
    dt[, Pathway := factor(Pathway, levels = c("SS", "SD", "DS", "DD"))]
    
    # Get Unit for Label
    unit <- attr(x$gen_intervals, "unit")
    unit_label <- switch(unit, 
                        "year" = "Years", "month" = "Months", "day" = "Days", 
                        "hour" = "Hours", "Years")
    
    # Basic formula
    f <- Mean ~ Pathway
    if ("Group" %in% names(dt)) {
      f <- Mean ~ Pathway | Group
    }
    
    p <- lattice::barchart(
      f, 
      data = dt,
      ylab = sprintf("Generation Interval (%s)", unit_label),
      xlab = "Pathway",
      main = "Generation Intervals by Pathway",
      col = "steelblue",
      border = "transparent",
      scales = list(x = list(rot = 0)),
      panel = function(x, y, ...) {
        lattice::panel.barchart(x, y, ...)
        lattice::panel.grid(h = -1, v = 0, col = "gray90", lty = "dotted")
      },
      ...
    )
    
    return(p)
    
  } else if (type == "ecg") {
    if (is.null(x$ecg)) {
      stop("No ECG data found in pedstats object.")
    }
    
    valid_metrics <- c("ECG", "FullGen", "MaxGen")
    if (!metric %in% valid_metrics) {
      stop("Invalid metric. Choose from: ", paste(valid_metrics, collapse = ", "))
    }
    
    # Check if metric column exists
    if (!metric %in% names(x$ecg)) {
      stop(sprintf("Metric '%s' not found in ECG data.", metric))
    }
    
    # Dynamic formula
    f <- as.formula(paste("~", metric))
    
    # Determine appropriate number of bins based on data range
    metric_data <- x$ecg[[metric]]
    n_unique <- length(unique(metric_data[!is.na(metric_data)]))
    n_bins <- min(30, max(10, ceiling(sqrt(length(metric_data)))))
    
    # For integer metrics (FullGen, MaxGen), use integer breaks
    if (metric %in% c("FullGen", "MaxGen")) {
      max_val <- max(metric_data, na.rm = TRUE)
      breaks <- seq(0, max_val + 1, by = 1)
      n_bins <- NULL
    } else {
      breaks <- NULL
    }
    
    p <- lattice::histogram(
      f, 
      data = x$ecg,
      xlab = metric,
      main = paste("Distribution of", metric),
      col = "lightgray",
      border = "white",
      breaks = breaks,
      nint = n_bins,
      type = "count",
      ...
    )
    
    return(p)
  }
}

#' @rdname vispstat
#' @export
plot.pedstats <- function(x, ...) {
  vispstat(x, ...)
}
