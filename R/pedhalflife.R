#' Calculate Information-Theoretic Diversity Half-Life
#'
#' @description
#' Calculates the diversity half-life (\eqn{T_{1/2}}) of a pedigree across
#' time points using a Renyi-2 entropy cascade framework.  The total loss rate of
#' genetic diversity is partitioned into three additive components:
#' \itemize{
#'   \item \eqn{\lambda_e}: foundation bottleneck (unequal founder
#'     contributions).
#'   \item \eqn{\lambda_b}: breeding bottleneck (overuse of key ancestors).
#'   \item \eqn{\lambda_d}: genetic drift (random sampling loss).
#' }
#'
#' The function rolls over time points defined by \code{timevar}, computing
#' \eqn{f_e} and \eqn{f_a} (via \code{\link{pedcontrib}}) and \eqn{f_g}
#' (via the internal coancestry engine) for each time point.  No redundant
#' Ne calculations are performed.
#'
#' @param ped A \code{tidyped} object.
#' @param timevar Character.
#'   Column name in \code{ped} that defines time points
#'   (e.g. \code{"Gen"}, \code{"Year"}).  Default: \code{"Gen"}.
#' @param at Optional vector of values selecting which time points to include
#'   (e.g., \code{2:4}, \code{2010:2020}, or non-consecutive \code{c(2015, 2018, 2022)}).
#'   Values must match entries in the \code{timevar} column.  Non-numeric values
#'   are accepted but the OLS time axis will fall back to sequential indices.
#'   If \code{NULL} (default), all non-NA unique values in \code{timevar} are used.
#' @param nsamples Integer.
#'   Sample size per time point for coancestry estimation
#'   (passed to the internal coancestry engine).  Default: 1000.
#' @param ncores Integer.
#'   Number of OpenMP threads for C++ backends.  Default: 1.
#' @param seed Integer or \code{NULL}.
#'   Random seed for reproducible coancestry sampling.
#'   Default: \code{NULL}.
#'
#' @return A list of class \code{pedhalflife} with two \code{data.table}
#'   components:
#'   \describe{
#'     \item{\code{timeseries}}{Per-time-point tracking with columns
#'       \code{Time} (time-point label from \code{timevar}), \code{NRef},
#'       \code{fe}, \code{fa}, \code{fg} and their log transformations
#'       (\code{lnfe}, \code{lnfa}, \code{lnfg}, \code{lnfafe},
#'       \code{lnfgfa}), plus \code{TimeStep} (numeric OLS time axis).}
#'     \item{\code{decay}}{Single-row table with \code{lambda_e},
#'       \code{lambda_b}, \code{lambda_d}, \code{lambda_total}, and
#'       \code{thalf}.}
#'   }
#'
#' @details
#' The mathematical identity underlying the cascade is:
#' \deqn{\ln f_g = \ln f_e + \ln(f_a / f_e) + \ln(f_g / f_a)}
#' Taking the negative time-slope of each term gives the \eqn{\lambda}
#' components, which sum exactly by linearity of OLS:
#' \deqn{\lambda_{total} = \lambda_e + \lambda_b + \lambda_d}
#'
#' \eqn{T_{1/2} = \ln 2 / \lambda_{total}} is the number of time-units
#' (time points, years, generations) for diversity to halve.
#'
#' @examples
#' \donttest{
#' library(visPedigree)
#' data(simple_ped)
#' tp <- tidyped(simple_ped)
#'
#' # 1. Calculate half-life using all available generations
#' hl <- pedhalflife(tp, timevar = "Gen")
#' print(hl)
#'
#' # 2. View the underlying log-linear decay plot
#' plot(hl, type = "log")
#'
#' # 3. Calculate half-life for a specific time window (e.g., Generations 2 to 4)
#' hl_subset <- pedhalflife(tp, timevar = "Gen", at = c(2, 3, 4))
#' print(hl_subset)
#' }
#'
#' @seealso \code{\link{pediv}}, \code{\link{pedcontrib}}, \code{\link{tidyped}}
#' @export
pedhalflife <- function(ped, timevar = "Gen", at = NULL, nsamples = 1000,
                        ncores = 1, seed = NULL) {
  ped <- ensure_complete_tidyped(ped, "pedhalflife()")

  if (!timevar %in% names(ped)) {
    stop(sprintf("Column '%s' not found in pedigree.", timevar))
  }

  # ---- Determine time points ----
  if (is.null(at)) {
    at <- sort(unique(ped[[timevar]]))
    at <- at[!is.na(at)]
  } else {
    at <- sort(unique(at))
  }

  if (length(at) < 2) {
    stop("At least two distinct time points are required to calculate half-life rates.")
  }

  # ---- One-time pre-computation: ECG on full pedigree ----
  ped_dt <- data.table::copy(ped)
  if (!"ECG" %in% names(ped_dt)) {
    ecg_dt <- pedecg(ped_dt)
    ped_dt <- merge(ped_dt, ecg_dt[, .(Ind, ECG)], by = "Ind", all.x = TRUE)
  }
  
  # Strip tidyped class from the working copy to avoid strict subsetting warnings
  data.table::setattr(ped_dt, "class", c("data.table", "data.frame"))

  # Configure OpenMP threads once (used by C++ coancestry backend)
  if (ncores > 1 && cpp_openmp_available()) {
    cpp_set_num_threads(ncores)
  }

  by_all <- ".CohortAll"
  ped_dt[[by_all]] <- "All"

  ts_list <- vector("list", length(at))

  message(sprintf("Calculating diversity across %d time points...", length(at)))
  for (i in seq_along(at)) {
    c_val <- at[i]
    ref_ids <- ped_dt[get(timevar) == c_val, Ind]

    if (length(ref_ids) == 0) next

    # ---- fe / fa via pedcontrib (public API, C++ backend) ----
    contrib <- suppressMessages(
      pedcontrib(ped, reference = ref_ids, mode = "both", top = 20)
    )
    s <- contrib$summary

    # ---- fg via internal coancestry engine (avoids redundant Ne calcs) ----
    fg_val <- NA_real_
    tryCatch({
      ped_subset <- as.data.table(ped_dt)[Ind %in% ref_ids]
      raw <- calc_ne_coancestry(ped_subset, ped_dt, by_all, nsamples,
                                seed = seed)
      if (!is.null(raw) && nrow(raw) > 0 && "fg" %in% names(raw)) {
        fg_val <- raw$fg[1]
      }
    }, error = function(e) {
      warning(sprintf("fg calculation failed for time point %s: %s",
                      as.character(c_val), e$message))
    })

    ts_list[[i]] <- data.table::data.table(
      Time   = c_val,
      NRef   = s$n_ref,
      fe     = s$f_e,
      fa     = s$f_a,
      fg     = fg_val
    )
  }

  ts_dt <- data.table::rbindlist(ts_list)

  if (nrow(ts_dt) < 2) {
    stop("Not enough valid time points to compute decay rates.")
  }

  # ---- Guard: require positive finite fe/fa/fg for log ----
  valid <- !is.na(ts_dt$fe) & !is.na(ts_dt$fa) & !is.na(ts_dt$fg) &
    ts_dt$fe > 0 & ts_dt$fa > 0 & ts_dt$fg > 0

  n_dropped <- sum(!valid)
  if (n_dropped > 0) {
    warning(sprintf(
      "%d time point(s) with NA or non-positive fe/fa/fg removed from decay analysis.",
      n_dropped
    ))
  }

  ts_valid <- ts_dt[valid]

  if (nrow(ts_valid) < 2) {
    stop("Not enough time points with valid positive fe/fa/fg to compute decay rates.")
  }

  # ---- Logarithmic metrics ----
  ts_valid[, lnfe := log(fe)]
  ts_valid[, lnfa := log(fa)]
  ts_valid[, lnfg := log(fg)]
  ts_valid[, lnfafe := log(fa / fe)]
  ts_valid[, lnfgfa := log(fg / fa)]

  # ---- Time axis ----
  t_vals <- suppressWarnings(as.numeric(ts_valid$Time))
  if (any(is.na(t_vals))) {
    warning("Time variable could not be parsed as numeric. Using sequential indices.")
    t_vals <- seq_len(nrow(ts_valid))
  }
  ts_valid[, TimeStep := t_vals]

  # ---- Decay rates (negative OLS slope of log vs time) ----
  calc_lambda <- function(y, t) {
    ok <- is.finite(y) & is.finite(t)
    if (sum(ok) < 2) return(NA_real_)
    y <- y[ok]; t <- t[ok]
    if (stats::sd(t) == 0) return(NA_real_)
    fit <- stats::lm(y ~ t)
    -stats::coef(fit)[[2]]
  }

  lambda_e     <- calc_lambda(ts_valid$lnfe,   ts_valid$TimeStep)
  lambda_b     <- calc_lambda(ts_valid$lnfafe, ts_valid$TimeStep)
  lambda_d     <- calc_lambda(ts_valid$lnfgfa, ts_valid$TimeStep)
  lambda_total <- calc_lambda(ts_valid$lnfg,   ts_valid$TimeStep)

  thalf <- if (!is.na(lambda_total) && lambda_total > 0) {
    log(2) / lambda_total
  } else {
    NA_real_
  }

  decay_dt <- data.table::data.table(
    LambdaE     = lambda_e,
    LambdaB     = lambda_b,
    LambdaD     = lambda_d,
    LambdaTotal = lambda_total,
    THalf       = thalf
  )

  result <- list(
    timeseries = ts_valid[],
    decay      = decay_dt[],
    timevar    = timevar
  )

  class(result) <- "pedhalflife"
  return(result)
}


# ---- S3 methods --------------------------------------------------------

#' @rdname pedhalflife
#' @param x A \code{pedhalflife} object.
#' @param ... Additional arguments (ignored).
#' @export
print.pedhalflife <- function(x, ...) {
  d <- x$decay

  fmt <- function(val) {
    if (is.na(val)) "        NA" else sprintf("%10.6f", val)
  }

  cat("Information-Theoretic Diversity Half-Life\n")
  cat("-----------------------------------------\n")
  cat(sprintf("Total Loss Rate (lambda_total): %s\n", fmt(d$LambdaTotal)))
  cat(sprintf("  Foundation  (lambda_e)      : %s\n", fmt(d$LambdaE)))
  cat(sprintf("  Bottleneck  (lambda_b)      : %s\n", fmt(d$LambdaB)))
  cat(sprintf("  Drift       (lambda_d)      : %s\n", fmt(d$LambdaD)))
  cat("-----------------------------------------\n")

  tv <- if (!is.null(x$timevar)) x$timevar else "time steps"
  if (!is.na(d$THalf) && d$THalf > 0) {
    cat(sprintf("Diversity Half-life (T_1/2)   : %10.2f (%s)\n", d$THalf, tv))
  } else {
    cat("Diversity Half-life (T_1/2)   : NA\n")
  }
  cat(sprintf("\nTimeseries: %d time points\n", nrow(x$timeseries)))
  invisible(x)
}


#' @rdname pedhalflife
#' @param type Character.
#'   \code{"log"} for log-transformed values;
#'   \code{"raw"} for \eqn{f_e}, \eqn{f_a}, \eqn{f_g}.
#' @export
plot.pedhalflife <- function(x, type = c("log", "raw"), ...) {
  type <- match.arg(type)

  if (!requireNamespace("lattice", quietly = TRUE)) {
    stop("Package 'lattice' is required for plotting.")
  }

  dt <- x$timeseries

  thalf <- x$decay$THalf
  tv    <- if (!is.null(x$timevar)) x$timevar else "Time"

  if (type == "log") {
    plot_data <- data.frame(
      Time   = rep(dt$TimeStep, 3),
      Value  = c(dt$lnfe, dt$lnfa, dt$lnfg),
      Metric = factor(rep(c("ln(fe)", "ln(fa)", "ln(fg)"), each = nrow(dt)),
                       levels = c("ln(fe)", "ln(fa)", "ln(fg)"))
    )

    # Pre-fit OLS on lnfg for the panel closure
    ols_fit <- stats::lm(lnfg ~ TimeStep, data = dt)

    lattice::xyplot(
      Value ~ Time,
      groups = Metric,
      data   = plot_data,
      type   = "b",
      auto.key = list(columns = 3),
      xlab   = tv,
      ylab   = "Log Diversity",
      main   = "Information-Theoretic Diversity Decay",
      panel  = function(...) {
        lattice::panel.xyplot(...)
        # OLS regression line for ln(fg) total decay
        lattice::panel.abline(ols_fit, lty = 2, col = "black", lwd = 1.5)
        # Half-life vertical reference
        if (is.finite(thalf)) {
          lattice::panel.abline(v = thalf, lty = 3, col = "grey40")
          lattice::panel.text(
            x = thalf, y = max(plot_data$Value, na.rm = TRUE),
            labels = sprintf("T1/2 = %.1f", thalf),
            pos = 4, cex = 0.8, col = "grey30"
          )
        }
      },
      ...
    )

  } else {
    plot_data <- data.frame(
      Time   = rep(dt$TimeStep, 3),
      Value  = c(dt$fe, dt$fa, dt$fg),
      Metric = factor(rep(c("fe", "fa", "fg"), each = nrow(dt)),
                       levels = c("fe", "fa", "fg"))
    )

    lattice::xyplot(
      Value ~ Time,
      groups = Metric,
      data   = plot_data,
      type   = "b",
      auto.key = list(columns = 3),
      xlab   = tv,
      ylab   = "Equivalent Numbers",
      main   = "Diversity Loss over Time",
      ...
    )
  }
}
