#' Calculate Generation Intervals
#'
#' Computes the generation intervals for the four gametic pathways:
#' Sire to Son (SS), Sire to Daughter (SD), Dam to Son (DS), and Dam to Daughter (DD).
#' The generation interval is defined as the age of the parents at the birth of their offspring.
#'
#' @param ped A \code{tidyped} object.
#' @param timevar Character. The name of the column containing the \strong{birth date}
#'   (or hatch date) of each individual. The column must be one of:
#'   \itemize{
#'     \item \code{Date} or \code{POSIXct} (recommended).
#'     \item A date string parseable by \code{as.POSIXct} (e.g., \code{"2020-03-15"}).
#'       Use \code{format} for non-ISO strings.
#'     \item A numeric year (e.g., \code{2020}). Automatically converted to
#'       \code{Date} (\code{"YYYY-07-01"}) with a message. For finer precision,
#'       convert to \code{Date} beforehand.
#'   }
#'   If \code{NULL}, auto-detects columns named \code{"BirthYear"}, \code{"Year"},
#'   \code{"BirthDate"}, or \code{"Date"}.
#' @param unit Character. Output time unit for the interval:
#'   \code{"year"} (default), \code{"month"}, \code{"day"}, or \code{"hour"}.
#' @param format Character. Optional format string for parsing \code{timevar}
#'   when it contains non-standard date strings (e.g., \code{"\%d/\%m/\%Y"}
#'   for \code{"15/03/2020"}).
#' @param cycle Numeric. Optional target (designed) length of one generation
#'   cycle expressed in \code{unit}s. When provided, an additional column
#'   \code{GenEquiv} is appended to the result, defined as:
#'   \deqn{GenEquiv_i = \frac{\bar{L}_i}{L_{cycle}}}
#'   where \eqn{\bar{L}_i} is the observed mean interval for pathway \eqn{i} and
#'   \eqn{L_{cycle}} is \code{cycle}. A value > 1 means the observed
#'   interval exceeds the target cycle (lower breeding efficiency).
#'   Example: for Pacific white shrimp with a 180-day target cycle, set
#'   \code{unit = "day", cycle = 180}.
#' @param by Character. Optional grouping column (e.g., \code{"Breed"}, \code{"Farm"}).
#'   If provided, intervals are calculated within each group.
#'
#' @return A \code{data.table} with columns:
#' \itemize{
#'   \item \code{Group}: Grouping level (if \code{by} is used).
#'   \item \code{Pathway}: One of \code{"SS"}, \code{"SD"}, \code{"DS"}, \code{"DD"},
#'     \code{"SO"}, \code{"DO"}, or \code{"Average"}.
#'     SS/SD/DS/DD require offspring sex; SO (Sire-Offspring) and DO (Dam-Offspring)
#'     are computed from all parent-offspring pairs regardless of offspring sex.
#'   \item \code{N}: Number of parent-offspring pairs used.
#'   \item \code{Mean}: Average generation interval in \code{unit}.
#'   \item \code{SD}: Standard deviation of the interval.
#'   \item \code{GenEquiv}: (Optional) Generation equivalents based on \code{cycle}.
#' }
#'
#' @details
#' Parent-offspring pairs with zero or negative intervals are excluded from
#' the calculation because they typically indicate data entry errors or
#' insufficient time resolution. If many zero intervals are expected (e.g.,
#' when using \code{unit = "year"} with annual spawners), consider using a
#' finer time unit such as \code{"month"} or \code{"day"}.
#'
#' Numeric year columns (e.g., \code{2020}) are automatically converted to
#' \code{Date} by appending \code{"-07-01"} (mid-year) as a reasonable default.
#' For more precise results, convert to \code{Date} before calling this function.
#'
#' @examples
#' \donttest{
#' # ---- Basic usage with package dataset (numeric Year auto-converted) ----
#' tped <- tidyped(big_family_size_ped)
#' gi <- pedgenint(tped, timevar = "Year")
#' gi
#'
#' # ---- Generation equivalents with cycle ----
#' gi2 <- pedgenint(tped, timevar = "Year", cycle = 2)
#' gi2
#' }
#'
#' @export
pedgenint <- function(ped, timevar = NULL, unit = c("year", "month", "day", "hour"),
                      format = NULL, cycle = NULL, by = NULL) {
  # Input validation
  ped <- ensure_complete_tidyped(ped, "pedgenint()")
  unit <- match.arg(unit)
  
  # Auto-detect timevar
  if (is.null(timevar)) {
    candidates <- c("BirthYear", "Year", "birth_year", "year", "BirthDate", "Date")
    match_col <- intersect(names(ped), candidates)
    if (length(match_col) > 0) {
      timevar <- match_col[1]
      message(sprintf("Using '%s' as time variable.", timevar))
    } else {
      stop("Please specify 'timevar' (e.g., 'BirthYear' or 'BirthDate').")
    }
  }
  
  if (!timevar %in% names(ped)) {
    stop(sprintf("Column '%s' not found in pedigree.", timevar))
  }
  
  # Check grouping column
  if (!is.null(by) && !by %in% names(ped)) {
    stop(sprintf("Grouping column '%s' not found.", by))
  }

  # Parse time to numeric axis
  numeric_time <- .parse_to_numeric_time(ped[[timevar]], unit = unit, format = format)
  
  # Prepare working data
  cols <- c("Ind", "Sire", "Dam", "Sex")
  if (!is.null(by)) cols <- c(cols, by)
  
  dt <- ped[, ..cols]
  dt[, Time := numeric_time]
  
  # Key: Ind -> Time
  time_map <- setNames(dt$Time, dt$Ind)
  
  calc_pathway <- function(parent_col, offspring_sex_label) {
    valid_offspring <- dt[!is.na(Sex) & Sex == offspring_sex_label]
    valid_pairs <- valid_offspring[!is.na(get(parent_col))]
    
    if (nrow(valid_pairs) == 0) return(NULL)
    
    p_ids <- valid_pairs[[parent_col]]
    p_times <- time_map[p_ids]
    o_times <- valid_pairs$Time
    
    intervals <- o_times - p_times
    # Filter valid intervals (ignore zero/negative which might be data errors)
    valid_idx <- !is.na(intervals) & intervals > 0
    
    if (sum(valid_idx) == 0) return(NULL)
    
    path_code <- paste0(
      ifelse(parent_col == "Sire", "S", "D"),
      ifelse(offspring_sex_label == "male", "S", "D")
    )
    
    res <- data.table(
      Pathway = path_code,
      Interval = intervals[valid_idx]
    )
    
    if (!is.null(by)) {
      res[, Group := valid_pairs[[by]][valid_idx]]
    }
    
    return(res)
  }
  
  res_list <- list(
    calc_pathway("Sire", "male"),   # SS
    calc_pathway("Sire", "female"), # SD
    calc_pathway("Dam", "male"),    # DS
    calc_pathway("Dam", "female")   # DD
  )
  
  all_intervals <- rbindlist(res_list)
  
  # Compute Average from ALL parent-offspring pairs regardless of offspring sex.
  # This avoids severely underestimating N when most offspring lack Sex info.
  calc_all_pathway <- function(parent_col) {
    valid_pairs <- dt[!is.na(get(parent_col))]
    if (nrow(valid_pairs) == 0) return(NULL)
    p_ids <- valid_pairs[[parent_col]]
    p_times <- time_map[p_ids]
    o_times <- valid_pairs$Time
    intervals <- o_times - p_times
    valid_idx <- !is.na(intervals) & intervals > 0
    if (sum(valid_idx) == 0) return(NULL)
    res <- data.table(
      Pathway = ifelse(parent_col == "Sire", "SO", "DO"),
      Interval = intervals[valid_idx]
    )
    if (!is.null(by)) {
      res[, Group := valid_pairs[[by]][valid_idx]]
    }
    return(res)
  }
  
  all_parent_intervals <- rbindlist(list(
    calc_all_pathway("Sire"),
    calc_all_pathway("Dam")
  ))
  
  if (nrow(all_intervals) == 0 && nrow(all_parent_intervals) == 0) {
    warning("No valid parent-offspring pairs found with birth dates.")
    res <- data.table(Pathway = character(0), N = integer(0), Mean = numeric(0), SD = numeric(0))
    attr(res, "unit") <- unit
    return(res)
  }
  
  agg_cols <- "Pathway"
  if (!is.null(by)) agg_cols <- c("Group", "Pathway")
  
  summary_dt <- if (nrow(all_intervals) > 0) {
    all_intervals[, .(
      N = .N,
      Mean = mean(Interval),
      SD = sd(Interval)
    ), by = agg_cols]
  } else {
    data.table(Pathway = character(0), N = integer(0), Mean = numeric(0), SD = numeric(0))
  }
  
  # SO/DO: Sire-Offspring and Dam-Offspring (sex-independent)
  so_do_agg_cols <- if (!is.null(by)) c("Group", "Pathway") else "Pathway"
  
  so_do_dt <- if (nrow(all_parent_intervals) > 0) {
    all_parent_intervals[, .(
      N = .N,
      Mean = mean(Interval),
      SD = sd(Interval)
    ), by = so_do_agg_cols]
  } else {
    data.table(Pathway = character(0), N = integer(0), Mean = numeric(0), SD = numeric(0))
  }
  
  # Average: all parent-offspring pairs combined
  overall_agg_cols <- if (!is.null(by)) "Group" else NULL
  
  overall_dt <- all_parent_intervals[, .(
    Pathway = "Average",
    N = .N,
    Mean = mean(Interval),
    SD = sd(Interval)
  ), by = overall_agg_cols]
  
  final_res <- rbind(summary_dt, so_do_dt, overall_dt, fill = TRUE)
  
  # Handle Cycle Length for Generation Equivalents
  if (!is.null(cycle) && is.numeric(cycle)) {
    final_res[, GenEquiv := Mean / cycle]
  }
  
  if (!is.null(by)) {
    setorderv(final_res, c("Group", "Pathway"))
  } else {
    setorder(final_res, Pathway)
  }
  
  attr(final_res, "unit") <- unit
  return(final_res[])
}

.parse_to_numeric_time <- function(x, unit = "year", format = NULL) {
  # Seconds per unit
  factors <- c(year = 365.25 * 86400, month = 30.4375 * 86400,
               day = 86400, hour = 3600)
  if (!unit %in% names(factors)) {
    stop("unit must be one of: ", paste(names(factors), collapse = ", "))
  }
  target_factor <- factors[unit]

  # --- Numeric input: assume integer years, auto-convert to Date ---
  if (is.numeric(x)) {
    message("Numeric time column detected. Converting to Date (YYYY-07-01). ",
            "For finer precision, convert to Date beforehand.")
    date_strings <- ifelse(is.na(x), NA_character_, paste0(as.integer(x), "-07-01"))
    x <- as.Date(date_strings)
  }

  # --- Date input ---
  if (inherits(x, "Date")) {
    seconds <- as.numeric(as.POSIXct(x, tz = "UTC"))
    return(seconds / target_factor)
  }

  # --- POSIXct input ---
  if (inherits(x, "POSIXct")) {
    return(as.numeric(x) / target_factor)
  }

  # --- Character input: parse date strings ---
  if (!is.character(x)) {
    stop("timevar must be Date, POSIXct, character, or numeric year. Got ",
         class(x)[1], ".")
  }

  u_x <- unique(x[!is.na(x)])
  if (length(u_x) == 0) return(rep(NA_real_, length(x)))

  parsed <- if (!is.null(format)) {
    as.POSIXct(u_x, format = format, tz = "UTC")
  } else {
    # Try ISO / common formats
    p <- suppressWarnings(as.POSIXct(u_x, tz = "UTC", optional = TRUE))
    if (all(is.na(p))) {
      p <- as.POSIXct(suppressWarnings(as.Date(u_x, optional = TRUE)), tz = "UTC")
    }
    if (all(is.na(p))) {
      p <- as.POSIXct(suppressWarnings(
        as.Date(u_x, format = "%Y/%m/%d", optional = TRUE)
      ), tz = "UTC")
    }
    # Fallback: bare year strings like "2020"
    if (all(is.na(p))) {
      num_try <- suppressWarnings(as.numeric(u_x))
      if (!all(is.na(num_try))) {
        p <- as.POSIXct(as.Date(paste0(as.integer(num_try), "-07-01")), tz = "UTC")
        message("Bare year strings detected. Converting to Date (YYYY-07-01). ",
                "For finer precision, convert to Date beforehand.")
      }
    }
    p
  }

  if (all(is.na(parsed))) {
    stop("Failed to parse time variable. Provide Date, POSIXct, or specify ",
         "'format' (e.g. '%Y-%m-%d').")
  }

  val_map <- setNames(as.numeric(parsed), u_x)
  numeric_seconds <- val_map[as.character(x)]

  return(numeric_seconds / target_factor)
}

#' Pedigree Subpopulations
#'
#' Summarizes pedigree subpopulations and group structure.
#'
#' When \code{by = NULL}, this function is a lightweight summary wrapper around
#' \code{\link{splitped}}, returning one row per disconnected pedigree component
#' plus an optional \code{"Isolated"} row for individuals with no known parents
#' and no offspring. When \code{by} is provided, it instead summarizes the pedigree
#' directly by the specified column (e.g. \code{"Gen"}, \code{"Year"}, \code{"Breed"}).
#'
#' @param ped A \code{tidyped} object.
#' @param by Character. The name of the column to group by. 
#'   If NULL, summarizes disconnected components via \code{\link{splitped}}.
#'
#' @return A \code{data.table} with columns:
#' \itemize{
#'   \item \code{Group}: Subpopulation label.
#'   \item \code{N}: Total individuals.
#'   \item \code{N_Sire}: Number of distinct sires.
#'   \item \code{N_Dam}: Number of distinct dams.
#'   \item \code{N_Founder}: Number of founders (parents unknown).
#' }
#'
#' @details
#' Use \code{pedsubpop()} when you want a compact analytical summary table.
#' Use \code{\link{splitped}} when you need the actual re-tidied sub-pedigree
#' objects for downstream plotting or analysis.
#'
#' @examples
#' tp <- tidyped(simple_ped)
#'
#' # Summarize disconnected pedigree components
#' pedsubpop(tp)
#'
#' # Summarize by an existing grouping variable
#' pedsubpop(tp, by = "Gen")
#' 
#' @export
pedsubpop <- function(ped, by = NULL) {
  ped <- ensure_tidyped(ped)
  
  if (is.null(by)) {
    splits <- splitped(ped)
    groups <- names(splits)
    n_grps <- length(splits)
    iso <- attr(splits, "isolated")
    n_iso <- length(iso)
    
    res_list <- vector("list", n_grps + (n_iso > 0))
    
    for (i in seq_along(splits)) {
      gped <- splits[[i]]
      res_list[[i]] <- list(
        Group = groups[i],
        N = nrow(gped),
        N_Sire = data.table::uniqueN(gped$Sire, na.rm = TRUE),
        N_Dam = data.table::uniqueN(gped$Dam, na.rm = TRUE),
        N_Founder = sum(is.na(gped$Sire) & is.na(gped$Dam))
      )
    }
    
    if (n_iso > 0) {
      res_list[[n_grps + 1]] <- list(
        Group = "Isolated",
        N = n_iso,
        N_Sire = 0,
        N_Dam = 0,
        N_Founder = n_iso
      )
    }
    
  } else {
    if (!by %in% names(ped)) stop(sprintf("Column '%s' not found.", by))
    
    res_list <- list(ped[, .(
      N = .N,
      N_Sire = data.table::uniqueN(Sire, na.rm = TRUE),
      N_Dam = data.table::uniqueN(Dam, na.rm = TRUE),
      N_Founder = sum(is.na(Sire) & is.na(Dam))
    ), by = by])
    
    setnames(res_list[[1]], by, "Group")
  }
  
  data.table::rbindlist(res_list)[]
}

#' Calculate Average Additive Genetic Relationship (\eqn{a_{ij}})
#'
#' Computes the average pairwise additive genetic relationship coefficients (\eqn{a_{ij}}) 
#' within cohorts or groups. The relationship \eqn{a_{ij}} is defined as twice the 
#' coancestry coefficient (\eqn{f_{ij}}), representing the expected proportion of 
#' genes shared by descent (e.g., 0.5 for full siblings).
#'
#' @param ped A \code{tidyped} object.
#' @param by Character. The column name to group by (e.g., "Year", "Breed", "Generation").
#' @param reference Character vector. An optional vector of reference individual IDs to calculate
#'   relationships for. If provided, only individuals matching these IDs in each group
#'   will be used. Default is NULL (use all individuals in the group).
#' @param compact Logical. Whether to use compact representation for large families to
#'   save memory. Recommended when pedigree size exceeds 25,000. Default is FALSE.
#'
#' @return A \code{data.table} with columns:
#' \itemize{
#'   \item A grouping identifier column, named after the \code{by} parameter (e.g., \code{Gen}, \code{Year}).
#'   \item \code{NTotal}: Total number of individuals in the group.
#'   \item \code{NUsed}: Number of individuals used in calculation (could be subset by reference).
#'   \item \code{MeanRel}: Average of off-diagonal elements in the Additive Relationship (A) matrix 
#'     for this group (\eqn{a_{ij} = 2f_{ij}}). Returns NA if the group has fewer than 2 individuals.
#' }
#' 
#' @examples
#' \donttest{
#' library(data.table)
#' # Use the sample dataset and simulate a birth year
#' tp <- tidyped(small_ped)
#' tp$Year <- 2010 + tp$Gen
#' 
#' # Example 1: Calculate average relationship grouped by Generation (default)
#' rel_by_gen <- pedrel(tp, by = "Gen")
#' print(rel_by_gen)
#' 
#' # Example 2: Calculate average relationship grouped by Year
#' rel_by_year <- pedrel(tp, by = "Year")
#' print(rel_by_year)
#' 
#' # Example 3: Filter calculations with a reference list in a chosen group
#' candidates <- c("N", "O", "P", "Q", "T", "U", "V", "X", "Y")
#' rel_subset <- pedrel(tp, by = "Gen", reference = candidates)
#' print(rel_subset)
#' }
#' 
#' @export
pedrel <- function(ped, by = "Gen", reference = NULL, compact = FALSE) {
  ped <- ensure_complete_tidyped(ped, "pedrel()")
  if (!by %in% names(ped)) stop(sprintf("Column '%s' not found.", by))
  
  groups <- unique(ped[[by]])
  groups <- groups[!is.na(groups)]
  
  res_list <- lapply(groups, function(g) {
    sub_ped_full <- ped[get(by) == g]
    n_total <- nrow(sub_ped_full)
    
    if (n_total < 2) {
      warning(sprintf("Group '%s' has less than 2 individuals, returning NA_real_.", g))
      return(data.table(Group = g, NTotal = n_total, NUsed = n_total, MeanRel = NA_real_))
    }
    
    if (!is.null(reference)) {
      sub_ped <- sub_ped_full[Ind %in% reference]
    } else {
      sub_ped <- sub_ped_full
    }
    n_used <- nrow(sub_ped)
    
    if (n_used < 2) {
      warning(sprintf("Group '%s' has less than 2 individuals after applying 'reference', returning NA_real_.", g))
      return(data.table(Group = g, NTotal = n_total, NUsed = n_used, MeanRel = NA_real_))
    }
    
    local_ped <- tryCatch({
      suppressMessages(tidyped(ped, cand = sub_ped$Ind, addnum = TRUE))
    }, error = function(e) {
      if (all(is.na(sub_ped$Sire) & is.na(sub_ped$Dam))) {
        return(NULL)
      }
      stop(e)
    })
    
    if (is.null(local_ped)) {
      mean_rel <- 0.0
    } else {
      target_inds <- sub_ped$Ind
      
      max_dense <- 25000L
      if (!compact && nrow(local_ped) > max_dense) {
        warning(sprintf(
          paste0("Group '%s': pedigree too large (%d individuals including ancestors) ",
                 "for dense A matrix (limit: %d). Returning NA.\n",
                 "Hint: use 'reference' parameter to reduce group size, or set compact = TRUE."),
          g, nrow(local_ped), max_dense))
        return(data.table(Group = g, NTotal = n_total, NUsed = n_used, MeanRel = NA_real_))
      }
      
      if (isTRUE(compact)) {
        A <- tryCatch({
          pedmat(local_ped, method = "A", sparse = FALSE, compact = TRUE)
        }, error = function(e) {
          warning(sprintf(
            paste0("Group '%s': %s\n",
                   "Hint: subset with 'reference' or use 'compact = TRUE'."),
            g, e$message))
          return(NULL)
        })
        
        if (is.null(A)) {
          return(data.table(Group = g, NTotal = n_total, NUsed = n_used, MeanRel = NA_real_))
        }
        
        c_map <- attr(A, "compact_map")
        c_map_target <- c_map[Ind %in% target_inds]
        
        freq_dt <- c_map_target[, .(W = .N), by = RepIndNum]
        W_idx <- freq_dt$RepIndNum
        W <- freq_dt$W
        
        A_rep <- A[W_idx, W_idx, drop = FALSE]
        sum_total_rep <- as.numeric(t(W) %*% (A_rep %*% W))
        sum_between <- sum_total_rep - sum(W^2 * diag(A_rep))
        
        sib_reps <- freq_dt[W > 1]
        sum_within <- 0.0
        
        if (nrow(sib_reps) > 0) {
          rep_parents <- unique(c_map_target[RepIndNum %in% sib_reps$RepIndNum,
                                             .(RepIndNum, rep_sire = SireNum, rep_dam = DamNum)])
          
          for (i in seq_len(nrow(sib_reps))) {
            w_i <- sib_reps$W[i]
            r_idx <- sib_reps$RepIndNum[i]
            sire_idx <- rep_parents[RepIndNum == r_idx, rep_sire]
            dam_idx <- rep_parents[RepIndNum == r_idx, rep_dam]
            A_sib <- 0.25 * (A[sire_idx, sire_idx] + A[dam_idx, dam_idx] + 2 * A[sire_idx, dam_idx])
            sum_within <- sum_within + w_i * (w_i - 1) * A_sib
          }
        }
        
        mean_rel <- (sum_between + sum_within) / (as.numeric(n_used) * (n_used - 1))
      } else {
        target_idx <- match(target_inds, local_ped$Ind)
        if (anyNA(target_idx)) {
          return(data.table(Group = g, NTotal = n_total, NUsed = n_used, MeanRel = NA_real_))
        }
        mean_rel <- tryCatch({
          cpp_mean_relationship(local_ped$SireNum, local_ped$DamNum, as.integer(target_idx))
        }, error = function(e) {
          warning(sprintf(
            paste0("Group '%s': %s\n",
                   "Hint: subset with 'reference' or use 'compact = TRUE'."),
            g, e$message))
          return(NA_real_)
        })
      }
    }
    
    data.table(Group = g, NTotal = n_total, NUsed = n_used, MeanRel = mean_rel)
  })
  
  result <- rbindlist(res_list)
  setnames(result, "Group", by)
  setorderv(result, cols = by)
  return(result[])
}

#' Calculate Equi-Generate Coefficient
#'
#' Estimates the number of distinct ancestral generations using the Equi-Generate Coefficient (ECG).
#' The ECG is calculated as 1/2 of the sum of the parents' ECG values plus 1.
#'
#' @param ped A \code{tidyped} object.
#'
#' @return A \code{data.table} with columns:
#' \itemize{
#'   \item \code{Ind}: Individual ID.
#'   \item \code{ECG}: Equi-Generate Coefficient.
#'   \item \code{FullGen}: Fully traced generations (depth of complete ancestry).
#'   \item \code{MaxGen}: Maximum ancestral generations (depth of deepest ancestor).
#' }
#'
#' @examples
#' tp <- tidyped(simple_ped)
#' ecg <- pedecg(tp)
#'
#' # ECG combines pedigree depth and completeness
#' head(ecg)
#'
#' # Individuals with deeper and more complete ancestry have larger ECG values
#' ecg[order(-ECG)][1:5]
#' 
#' @references
#' Boichard, D., Maignel, L., & Verrier, E. (1997). The value of using probabilities of gene origin to measure genetic variability in a population. Genetics Selection Evolution, 29(1), 5.
#'
#' @export
pedecg <- function(ped) {
  ped <- ensure_complete_tidyped(ped, "pedecg()")
  
  n <- nrow(ped)
  ecg <- numeric(n)
  full_gen <- numeric(n)
  max_gen_val <- numeric(n)
  
  sire_nums <- ped$SireNum
  dam_nums <- ped$DamNum
  
  gen_indices <- split(ped$IndNum, ped$Gen)
  gen_levels <- sort(as.integer(names(gen_indices)))
  
  for (g in gen_levels) {
    idx <- gen_indices[[as.character(g)]]
    k <- length(idx)
    
    sires <- sire_nums[idx]
    dams <- dam_nums[idx]
    
    sire_ecg <- numeric(k) 
    dam_ecg <- numeric(k)
    sire_full <- rep(-1, k)
    dam_full <- rep(-1, k)
    sire_max <- rep(-1, k)
    dam_max <- rep(-1, k)
    
    valid_s <- !is.na(sires) & sires > 0
    if (any(valid_s)) {
      parents <- sires[valid_s]
      sire_ecg[valid_s] <- ecg[parents]
      sire_full[valid_s] <- full_gen[parents]
      sire_max[valid_s] <- max_gen_val[parents]
    }
    
    valid_d <- !is.na(dams) & dams > 0
    if (any(valid_d)) {
      parents <- dams[valid_d]
      dam_ecg[valid_d] <- ecg[parents]
      dam_full[valid_d] <- full_gen[parents]
      dam_max[valid_d] <- max_gen_val[parents]
    }
    
    term_s <- ifelse(valid_s, sire_ecg + 1, 0)
    term_d <- ifelse(valid_d, dam_ecg + 1, 0)
    ecg[idx] <- 0.5 * (term_s + term_d)
    
    full_gen[idx] <- pmin(sire_full, dam_full) + 1
    max_gen_val[idx] <- pmax(sire_max, dam_max) + 1
  }
  
  res <- data.table(
    Ind = ped$Ind,
    ECG = ecg,
    FullGen = full_gen,
    MaxGen = max_gen_val
  )

  return(res[])
}

#' Pedigree Statistics
#'
#' Calculates comprehensive statistics for a pedigree, including population structure,
#' generation intervals, and ancestral depth.
#'
#' @param ped A \code{tidyped} object.
#' @param timevar Optional character. Name of the column containing the
#'   \strong{birth date} (or hatch date) of each individual.
#'   Accepted column formats:
#'   \itemize{
#'     \item \code{Date} or \code{POSIXct} (recommended).
#'     \item A date string parseable by \code{as.POSIXct}
#'       (e.g., \code{"2020-06-15"}).  Use \code{format} via \code{...}
#'       for non-ISO strings.
#'     \item A numeric year (e.g., \code{2020}).  Automatically converted
#'       to \code{Date} (\code{"YYYY-07-01"}) with a message.
#'   }
#'   If \code{NULL}, attempts auto-detection from common column names
#'   (\code{"BirthYear"}, \code{"Year"}, \code{"BirthDate"}, etc.).
#' @param unit Character. Time unit for reporting generation intervals:
#'   \code{"year"} (default), \code{"month"}, \code{"day"}, or \code{"hour"}.
#' @param cycle Numeric. Optional target generation cycle length in
#'   \code{unit}s. When provided, \code{gen_intervals} will include a
#'   \code{GenEquiv} column (observed Mean / cycle). See
#'   \code{\link{pedgenint}} for details.
#' @param ecg Logical. Whether to compute equivalent complete generations
#'   for each individual via \code{\link{pedecg}}. Default \code{TRUE}.
#' @param genint Logical. Whether to compute generation intervals via
#'   \code{\link{pedgenint}}. Requires a detectable \code{timevar} column.
#'   Default \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link{pedgenint}},
#'   e.g., \code{format} for custom date parsing or \code{by} for grouping.
#'
#' @return An object of class \code{pedstats}, which is a list containing:
#' \itemize{
#'   \item \code{summary}: A \code{data.table} with one row summarising the
#'     whole pedigree.  Columns:
#'     \itemize{
#'       \item \code{N} — total number of individuals.
#'       \item \code{NSire} — number of unique sires.
#'       \item \code{NDam} — number of unique dams.
#'       \item \code{NFounder} — number of founder individuals
#'         (both parents unknown).
#'       \item \code{MaxGen} — maximum generation number.
#'     }
#'   \item \code{ecg}: A \code{data.table} with one row per individual
#'     (\code{NULL} if \code{ecg = FALSE}).  Columns:
#'     \itemize{
#'       \item \code{Ind} — individual identifier.
#'       \item \code{ECG} — equivalent complete generations.
#'       \item \code{FullGen} — number of fully known generations.
#'       \item \code{MaxGen} — maximum traceable generation depth.
#'     }
#'   \item \code{gen_intervals}: A \code{data.table} of generation intervals
#'     (\code{NULL} if no \code{timevar} is detected or
#'     \code{genint = FALSE}).  Columns:
#'     \itemize{
#'       \item \code{Pathway} — gametic pathway label.  Seven values:
#'         \code{"SS"} (sire to son), \code{"SD"} (sire to daughter),
#'         \code{"DS"} (dam to son), \code{"DD"} (dam to daughter) —
#'         require offspring sex; \code{"SO"} (sire to offspring) and
#'         \code{"DO"} (dam to offspring) — sex-independent; and
#'         \code{"Average"} — all parent-offspring pairs combined.
#'       \item \code{N} — number of parent-offspring pairs.
#'       \item \code{Mean} — mean generation interval.
#'       \item \code{SD} — standard deviation of the interval.
#'       \item \code{GenEquiv} — \code{Mean / cycle} (only present when
#'         \code{cycle} is supplied).
#'     }
#' }
#'
#' @examples
#' \donttest{
#' # ---- Without time variable ----
#' tp <- tidyped(simple_ped)
#' ps <- pedstats(tp)
#' ps$summary
#' ps$ecg
#'
#' # ---- With annual Year column (big_family_size_ped) ----
#' tp2 <- tidyped(big_family_size_ped)
#' ps2 <- pedstats(tp2, timevar = "Year")
#' ps2$summary
#' ps2$gen_intervals
#' }
#'
#' @export
pedstats <- function(ped, timevar = NULL, unit = "year", cycle = NULL, ecg = TRUE, genint = TRUE, ...) {
  ped <- if (isTRUE(ecg) || isTRUE(genint)) {
    ensure_complete_tidyped(ped, "pedstats()")
  } else {
    ensure_tidyped(ped)
  }
  
  summ <- data.table(
    N = nrow(ped),
    NSire = data.table::uniqueN(ped$Sire, na.rm = TRUE),
    NDam = data.table::uniqueN(ped$Dam, na.rm = TRUE),
    NFounder = sum(is.na(ped$Sire) & is.na(ped$Dam)),
    MaxGen = max(ped$Gen)
  )
  
  ecg_dt <- NULL
  if (ecg) {
    ecg_dt <- pedecg(ped)
  }
  
  gen_int <- NULL
  
  target_timevar <- timevar
  if (is.null(target_timevar)) {
    candidates <- c("BirthYear", "Year", "birth_year", "year", "BirthDate", "Date")
    match_col <- intersect(names(ped), candidates)
    if (length(match_col) > 0) target_timevar <- match_col[1]
  }
  
  if (genint && !is.null(target_timevar) && target_timevar %in% names(ped)) {
    gen_int <- tryCatch({
      pedgenint(ped, timevar = target_timevar, unit = unit, cycle = cycle, ...)
    }, error = function(e) {
      warning("Failed to calculate generation intervals: ", e$message)
      NULL
    })
  }
  
  res <- list(
    summary = summ,
    ecg = ecg_dt,
    gen_intervals = gen_int
  )
  class(res) <- "pedstats"
  res
}

#' Print Pedigree Statistics
#'
#' @param x A \code{pedstats} object.
#' @param ... Additional arguments.
#'
#' @export
print.pedstats <- function(x, ...) {
  cat("Pedigree Statistics\n")
  cat("===================\n")
  
  s <- x$summary
  cat(sprintf("Total Individuals: %d\n", s$N))
  cat(sprintf("Sires: %d | Dams: %d | Founders: %d\n", s$NSire, s$NDam, s$NFounder))
  cat(sprintf("Max Generation: %d\n", s$MaxGen))
  
  if (!is.null(x$gen_intervals)) {
    unit <- attr(x$gen_intervals, "unit")
    unit_label <- switch(unit, 
                        "year" = "Years", "month" = "Months", "day" = "Days", 
                        "hour" = "Hours", "gen" = "Units", "Years")
    cat(sprintf("\nGeneration Intervals (%s):\n", unit_label))
    print(x$gen_intervals)
  }
  
  if (!is.null(x$ecg)) {
    cat("\nAncestral Depth (ECG) Summary:\n")
    print(summary(x$ecg$ECG))
  }
  
  invisible(x)
}

#' Calculate Effective Population Size
#'
#' Calculates the effective population size (Ne) based on the rate of
#' coancestry, the rate of inbreeding, or demographic parent numbers.
#'
#' @param ped A \code{tidyped} object.
#' @param method Character. The method to compute Ne. One of \code{"coancestry"} (default), \code{"inbreeding"}, or \code{"demographic"}.
#' @param by Character. The name of the column used to group cohorts (e.g., "Year", "BirthYear").
#'   If NULL, calculates overall Ne for all individuals.
#' @param reference Character vector. Optional subset of individual IDs defining the reference cohort.
#'   If NULL, uses all individuals in the pedigree.
#' @param nsamples Integer. Number of individuals to randomly sample per cohort when using the \code{"coancestry"} method. Very large cohorts will be sampled down to this size to save memory and time (default: 1000).
#' @param ncores Integer. Number of cores for parallel processing. Currently only effective for \code{method = "coancestry"} (default: 1).
#' @param seed Integer or NULL. Random seed passed to \code{set.seed()} before sampling
#'   in the coancestry method, ensuring reproducible \eqn{N_e} estimates.
#'   Default is \code{NULL} (no fixed seed).
#'
#' @return A \code{data.table} with columns:
#' \itemize{
#'   \item \code{Cohort}: Cohort or grouping variable value.
#'   \item \code{N}: Number of individuals in the cohort.
#'   \item \code{Ne}: Calculated effective population size.
#'   \item \code{...}: Additional columns depending on the selected method (e.g., \code{NSampled}, \code{DeltaC}, \code{MeanF}, \code{DeltaF}, \code{Nm}, \code{Nf}).
#' }
#' 
#' @details
#' The effective population size can be calculated using one of three methods:
#' 
#' \itemize{
#'   \item \strong{"coancestry"} (Default): Based on the rate of coancestry (\eqn{f_{ij}}) between 
#'   pairs of individuals. In this context, \eqn{f_{ij}} is the probability that two alleles 
#'   randomly sampled from individuals $i$ and $j$ are identical by descent, which is equivalent 
#'   to half the additive genetic relationship (\eqn{f_{ij} = a_{ij} / 2}). 
#'   This method is generally more robust as it accounts for full genetic drift and 
#'   bottlenecks (Cervantes et al., 2011).
#'   \deqn{\Delta c_{ij} = 1 - (1 - c_{ij})^{1/(\frac{ECG_i + ECG_j}{2})}}
#'   \deqn{N_e = \frac{1}{2 \overline{\Delta c}}}
#'   To handle large populations, this method samples \code{nsamples} individuals per cohort 
#'   and computes the mean rate of coancestry among them.
#'   
#'   \item \strong{"inbreeding"}: Based on the individual rate of inbreeding ($F_i$) (Gutiérrez et al., 2008, 2009).
#'   \deqn{\Delta F_i = 1 - (1 - F_i)^{1/(ECG_i - 1)}}
#'   \deqn{N_e = \frac{1}{2 \overline{\Delta F}}}
#'   
#'   \item \strong{"demographic"}: Based on the demographic census of breeding males and females (Wright, 1931).
#'   \deqn{N_e = \frac{4 N_m N_f}{N_m + N_f}}
#'   Where \eqn{N_m} and \eqn{N_f} are the number of unique male and female parents contributing to the cohort.
#' }
#' 
#' @references
#' Cervantes, I., Goyache, F., Molina, A., Valera, M., & Gutiérrez, J. P. (2011). Estimation of effective population size from the rate of coancestry in pedigreed populations. \emph{Journal of Animal Breeding and Genetics}, 128(1), 56-63.
#' 
#' Gutiérrez, J. P., Cervantes, I., Molina, A., Valera, M., & Goyache, F. (2008). Individual increase in inbreeding allows estimating effective sizes from pedigrees. \emph{Genetics Selection Evolution}, 40(4), 359-370.
#' 
#' Gutiérrez, J. P., Cervantes, I., & Goyache, F. (2009). Improving the estimation of realized effective population sizes in farm animals. \emph{Journal of Animal Breeding and Genetics}, 126(4), 327-332.
#' 
#' Wright, S. (1931). Evolution in Mendelian populations. \emph{Genetics}, 16(2), 97-159.
#'
#' @examples
#' \donttest{
#' # Coancestry-based Ne (default) using a simple pedigree grouped by year
#' tp_simple <- tidyped(simple_ped)
#' tp_simple$BirthYear <- 2000 + tp_simple$Gen
#' ne_coan <- suppressMessages(pedne(tp_simple, by = "BirthYear", seed = 42L))
#' ne_coan
#'
#' # Inbreeding-based Ne using an inbred pedigree
#' tp_inbred <- tidyped(inbred_ped)
#' ne_inb <- suppressMessages(pedne(tp_inbred, method = "inbreeding", by = "Gen"))
#' ne_inb
#'
#' # Demographic Ne from the number of contributing sires and dams
#' ne_demo <- suppressMessages(pedne(tp_simple, method = "demographic", by = "BirthYear"))
#' ne_demo
#' }
#' 
#' @export
pedne <- function(ped, method = c("coancestry", "inbreeding", "demographic"),
                  by = NULL, reference = NULL, nsamples = 1000, ncores = 1,
                  seed = NULL) {
  
  ped <- if (match.arg(method) %in% c("coancestry", "inbreeding")) {
    ensure_complete_tidyped(ped, sprintf("pedne(method = \"%s\")", match.arg(method)))
  } else {
    ensure_tidyped(ped)
  }
  method <- match.arg(method)
  
  # Protect the original data from being modified by reference
  ped_dt <- data.table::copy(ped)
  
  # Handle 'by' parameter
  if (is.null(by)) {
    message("No 'by' parameter specified. Calculating Ne for the entire dataset. Grouping by a cohort or year variable is often more biologically meaningful.")
    by <- ".CohortAll"
    ped_dt[[by]] <- "All"
  } else {
    if (!by %in% names(ped_dt)) {
      stop(sprintf("Column '%s' not found in pedigree.", by))
    }
  }
  
  # Configure threads
  if (ncores > 1) {
    if (cpp_openmp_available()) {
      cpp_set_num_threads(ncores)
    } else {
      warning("OpenMP not available. Using single core.")
    }
  }

  # Pre-calculate mandatory genetic parameters on FULL pedigree before subsetting
  if (method %in% c("inbreeding", "coancestry")) {
    if (method == "inbreeding" && !"f" %in% names(ped_dt)) {
      ped_dt <- inbreed(ped_dt)
    }
    if (!"ECG" %in% names(ped_dt)) {
      ecg_dt <- pedecg(ped_dt)
      ped_dt <- merge(ped_dt, ecg_dt[, .(Ind, ECG)], by = "Ind", all.x = TRUE)
    }
  }
  
  # Filter reference population if specified
  if (!is.null(reference)) {
    if (!all(reference %in% ped_dt$Ind)) {
      missing <- reference[!reference %in% ped_dt$Ind]
      warning(sprintf("Reference IDs not found in pedigree: %s", paste(missing, collapse = ", ")))
    }
    ped_subset <- ped_dt[Ind %in% reference]
  } else {
    ped_subset <- ped_dt
  }
  
  # Dispatch to specific calculator
  res <- switch(method,
    "inbreeding"  = calc_ne_inbreeding(ped_subset, by),
    "coancestry"  = {
      raw <- calc_ne_coancestry(ped_subset, ped_dt, by, nsamples, seed = seed)
      # Strip fg/MeanCoan/NSampledCoan: exposed only via pediv(), not pedne()
      cols_to_drop <- intersect(names(raw), c("fg", "MeanCoan", "NSampledCoan"))
      if (length(cols_to_drop) > 0) raw[, (cols_to_drop) := NULL]
      raw
    },
    "demographic" = calc_ne_demographic(ped_subset, by)
  )

  return(res[])
}

calc_ne_inbreeding <- function(ped, by) {
  # Inbreeding coefficients and ECG should already be calculated on the full pedigree
  if (!"f" %in% names(ped) || !"ECG" %in% names(ped)) {
    stop("f or ECG is missing in subset tracking.")
  }
  
  ped_work <- ped[ECG >= 1]
  
  if (nrow(ped_work) == 0) {
    warning("No individuals with sufficient pedigree depth (ECG >= 1).")
    return(data.table(Cohort = numeric(0), N = integer(0), Ne = numeric(0)))
  }

  ped_work[, CohortLabel := get(by)]
  ped_work[, DeltaF_Ind := ifelse(ECG > 1, 1 - (1 - f)^(1 / (ECG - 1)), NA_real_)]
  
  # Filter valid DeltaF
  ped_valid <- ped_work[!is.na(DeltaF_Ind) & is.finite(DeltaF_Ind) & DeltaF_Ind >= 0]
  
  if (nrow(ped_valid) == 0) {
    warning("No valid DeltaF values calculated.")
    return(data.table(Cohort = numeric(0), N = integer(0), Ne = numeric(0)))
  }

  result <- ped_valid[, .(
    N = .N,
    MeanF = mean(f, na.rm = TRUE),
    DeltaF = mean(DeltaF_Ind, na.rm = TRUE),
    MeanECG = mean(ECG, na.rm = TRUE)
  ), by = .(Cohort = CohortLabel)]
  
  result[, Ne := ifelse(DeltaF > 0, 1 / (2 * DeltaF), NA_real_)]
  setorder(result, Cohort)
  
  return(result[])
}

calc_ne_demographic <- function(ped, by) {
  # Ensure numbers are present
  if (!all(c("Sire", "Dam", "Sex") %in% names(ped))) {
     stop("Pedigree must contain 'Sire', 'Dam', and 'Sex' columns for demographic Ne.")
  }

  ped <- copy(ped)
  ped[, CohortLabel := get(by)]
  
  # Need to find parents contributing to EACH cohort
  # Parent IDs are in Sire and Dam columns
  # We should count how many UNIQUE Sires and Dams produced offspring in this cohort
  
  result <- ped[, {
    sires <- unique(Sire[!is.na(Sire) & Sire != "" & Sire != "0"])
    dams <- unique(Dam[!is.na(Dam) & Dam != "" & Dam != "0"])
    
    nm <- length(sires)
    nf <- length(dams)
    
    ne_val <- if (nm > 0 && nf > 0) (4 * nm * nf) / (nm + nf) else NA_real_
    
    list(
      N = .N,
      Nm = nm,
      Nf = nf,
      Ne = ne_val
    )
  }, by = .(Cohort = CohortLabel)]
  
  setorder(result, Cohort)
  return(result[])
}

calc_ne_coancestry <- function(ped_subset, ped_full, by, nsamples,
                               seed = NULL, max_trace = 25000L) {
  # ECG is already merged into ped_subset from ped_full by the caller function.
  # Returns extended table with fg/MeanCoan/NSampledCoan columns.
  # NOTE: pedne() strips fg columns before exposing its public output;
  #       pediv() calls this function directly to access fg without duplication.

  ped_subset[, CohortLabel := get(by)]
  cohorts <- sort(unique(ped_subset$CohortLabel))

  results_list <- list()

  for (coh in cohorts) {
    sub_ped <- ped_subset[CohortLabel == coh]
    n_total <- nrow(sub_ped)

    if (n_total == 0) next

    # Reproducible sampling
    if (!is.null(seed)) set.seed(seed)

    # Initial sampling
    if (n_total > nsamples) {
      sampled_inds <- sample(sub_ped$Ind, nsamples)
    } else {
      sampled_inds <- sub_ped$Ind
    }

    # Trace back ancestors using the FULL pedigree
    trace_ped <- tidyped(ped_full, cand = sampled_inds, trace = "up",
                         addgen = TRUE, addnum = TRUE)

    # Adaptive resampling: if traced pedigree exceeds memory limit, reduce sample
    if (nrow(trace_ped) > max_trace) {
      n_reduced <- max(2L, floor(length(sampled_inds) * max_trace / nrow(trace_ped)))
      if (!is.null(seed)) set.seed(seed)
      sampled_inds <- sample(sub_ped$Ind, min(n_reduced, n_total))
      trace_ped <- tidyped(ped_full, cand = sampled_inds, trace = "up",
                           addgen = TRUE, addnum = TRUE)
    }

    # Must be sorted by IndNum for correct C++ indexing (1..N)
    setorder(trace_ped, IndNum)

    # Re-identify target numeric IDs (1-based)
    target_nums <- trace_ped[Ind %in% sampled_inds, IndNum]

    # Ensure ECG is present (propagate from full pedigree if missing after trace)
    if (!"ECG" %in% names(trace_ped)) {
      ecg_map <- ped_full[Ind %in% trace_ped$Ind, .(Ind, ECG)]
      trace_ped <- merge(trace_ped, ecg_map, by = "Ind", all.x = TRUE)
      setorder(trace_ped, IndNum)
    }
    trace_ped[is.na(ECG), ECG := 0]

    # ---- 1. DeltaC / Ne ----
    avg_delta_c <- tryCatch({
      cpp_calculate_sampled_coancestry_delta(
        trace_ped$SireNum,
        trace_ped$DamNum,
        target_nums,
        trace_ped$ECG
      )
    }, error = function(e) {
      warning(paste("Error in coancestry calculation for cohort", coh, ":", e$message))
      NA_real_
    })

    ne_val <- if (!is.na(avg_delta_c) && avg_delta_c > 0) 1 / (2 * avg_delta_c) else NA_real_

    # ---- 2. fg: founder genome equivalents (Caballero & Toro 2000) ----
    # Reuses the same traced pedigree — zero extra tidyped() cost.
    # Corrected mean coancestry:
    #   C_bar = (N-1)/N * rbar_off/2 + (1 + mean_F_s) / (2N)
    # where N = full reference cohort size (NOT sample size).
    fg_val     <- NA_real_
    mean_coan  <- NA_real_

    tryCatch({
      # Off-diagonal mean relationship (cpp returns off-diagonal only)
      rbar_off <- cpp_mean_relationship(
        trace_ped$SireNum,
        trace_ped$DamNum,
        target_nums
      )

      # Inbreeding for diagonal correction (on traced pedigree — fast)
      f_traced  <- cpp_calculate_inbreeding(trace_ped$SireNum, trace_ped$DamNum)$f
      mean_f_s  <- mean(f_traced[target_nums], na.rm = TRUE)

      # Diagonal-corrected mean coancestry weighted by full N
      N_ref     <- n_total
      mean_coan <- (N_ref - 1) / N_ref * rbar_off / 2 +
                   (1 + mean_f_s) / (2 * N_ref)

      fg_val    <- if (!is.na(mean_coan) && mean_coan > 0) 1 / (2 * mean_coan) else NA_real_
    }, error = function(e) {
      warning(paste("Error in fg calculation for cohort", coh, ":", e$message))
    })

    results_list[[as.character(coh)]] <- data.table(
      Cohort       = coh,
      N            = n_total,
      NSampled     = length(sampled_inds),
      DeltaC       = avg_delta_c,
      Ne           = ne_val,
      MeanCoan     = mean_coan,
      fg           = fg_val,
      NSampledCoan = length(sampled_inds)
    )
  }

  result <- rbindlist(results_list)
  return(result[])
}

# Legacy wrapper/stub for documentation or compatibility if needed
# (none required as we replaced pedne in place) -> THIS LINE INTENTIONALLY LEFT BLANK
fill_phantoms <- function(ped) {
  missing_sire <- is.na(ped$Sire) & !is.na(ped$Dam)
  missing_dam <- !is.na(ped$Sire) & is.na(ped$Dam)
  
  if (!any(missing_sire) && !any(missing_dam)) {
    return(ped)
  }
  
  ped_new <- data.table::copy(ped)
  new_rows <- list()
  
  if (any(missing_sire)) {
    phantoms_sire <- paste0("*Phantom_Sire_", ped_new$Ind[missing_sire])
    ped_new$Sire[missing_sire] <- phantoms_sire
    new_rows[[1]] <- data.table::data.table(
      Ind = phantoms_sire,
      Sire = NA_character_,
      Dam = NA_character_
    )
  }
  
  if (any(missing_dam)) {
    phantoms_dam <- paste0("*Phantom_Dam_", ped_new$Ind[missing_dam])
    ped_new$Dam[missing_dam] <- phantoms_dam
    new_rows[[2]] <- data.table::data.table(
      Ind = phantoms_dam,
      Sire = NA_character_,
      Dam = NA_character_
    )
  }
  
  new_dt <- data.table::rbindlist(new_rows, fill = TRUE)
  extra_cols <- setdiff(names(ped_new), names(new_dt))
  if (length(extra_cols) > 0) {
    for (col in extra_cols) {
      data.table::set(new_dt, j = col, value = NA)
    }
  }
  
  ped_combined <- rbind(ped_new, new_dt, fill = TRUE)
  
  old_cols <- intersect(names(ped_combined), c("IndNum", "SireNum", "DamNum"))
  if (length(old_cols) > 0) {
    ped_combined[, (old_cols) := NULL]
  }
  
  return(tidyped(ped_combined, addnum = TRUE))
}

# Internal: Shannon entropy -> effective number (q=1 Hill number)
# p must be a probability vector (sums to 1); zeros are filtered out.
calc_shannon_effective <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0L) return(NA_real_)
  exp(-sum(p * log(p)))
}

#' Calculate Founder and Ancestor Contributions
#'
#' Calculates genetic contributions from founders and influential ancestors.
#' Implements the gene dropping algorithm for founder contributions and Boichard's 
#' algorithm for ancestor contributions to estimate the effective number of founders ($f_e$) 
#' and ancestors ($f_a$).
#'
#' @param ped A \code{tidyped} object.
#' @param reference Character vector. Optional subset of individual IDs defining the reference population.
#'   If NULL, uses all individuals in the most recent generation.
#' @param mode Character. Type of contribution to calculate:
#' \itemize{
#'   \item \code{"founder"}: Founder contributions ($f_e$).
#'   \item \code{"ancestor"}: Ancestor contributions ($f_a$).
#'   \item \code{"both"}: Both founder and ancestor contributions.
#' }
#' @param top Integer. Number of top contributors to return. Default is 20.
#'
#' @return A list with class \code{pedcontrib} containing:
#' \itemize{
#'   \item \code{founders}: A \code{data.table} of founder contributions (if mode includes "founder", or "both").
#'   \item \code{ancestors}: A \code{data.table} of ancestor contributions (if mode includes "ancestor", or "both").
#'   \item \code{summary}: A \code{list} of summary statistics including:
#'     \itemize{
#'       \item \code{f_e}: Classical effective number of founders (\eqn{q=2}, Lacy 1989).
#'       \item \code{f_e_H}: Information-theoretic effective number of founders
#'         (\eqn{q=1}, Shannon entropy): \eqn{f_e^{(H)} = \exp(-\sum p_i \ln p_i)}.
#'       \item \code{f_a}: Classical effective number of ancestors (\eqn{q=2}, Boichard 1997).
#'       \item \code{f_a_H}: Information-theoretic effective number of ancestors
#'         (\eqn{q=1}): \eqn{f_a^{(H)} = \exp(-\sum q_k \ln q_k)}.
#'     }
#' }
#' 
#' Each contribution table contains:
#' \itemize{
#'   \item \code{Ind}: Individual ID.
#'   \item \code{Contrib}: Contribution to the reference population (0-1).
#'   \item \code{CumContrib}: Cumulative contribution.
#'   \item \code{Rank}: Rank by contribution.
#' }
#' 
#' @details
#' **Founder Contributions ($f_e$):**
#' Calculated by probabilistic gene flow from founders to the reference cohort. 
#' When individual ancestors with one unknown parent exist, "phantom" parents are temporarily injected correctly conserving the probability mass.
#' 
#' **Ancestor Contributions ($f_a$):**
#' Calculated using Boichard's iterative algorithm (1997), accounting for:
#' \itemize{
#'   \item Marginal genetic contribution of each ancestor
#'   \item Long-term contributions through multiple pathways
#' }
#' The parameter $f_a$ acts as a stringent metric since it identifies the bottlenecks of genetic variation in pedigrees. 
#' 
#' @examples
#' \donttest{
#' library(data.table)
#' # Load a sample pedigree
#' tp <- tidyped(small_ped)
#'
#' # Calculate both founder and ancestor contributions for reference population
#' ref_ids <- c("Z1", "Z2", "X", "Y")
#' contrib <- pedcontrib(tp, reference = ref_ids, mode = "both")
#'
#' # Print results including f_e, f_e(H), f_a, and f_a(H)
#' print(contrib)
#'
#' # Access Shannon-entropy effective numbers directly
#' contrib$summary$f_e_H   # Information-theoretic effective founders (q=1)
#' contrib$summary$f_e     # Classical effective founders (q=2)
#' contrib$summary$f_a_H   # Information-theoretic effective ancestors (q=1)
#' contrib$summary$f_a     # Classical effective ancestors (q=2)
#'
#' # Diversity ratio rho > 1 indicates long-tail founder value
#' contrib$summary$f_e_H / contrib$summary$f_e
#' }
#'
#' @references
#' Boichard, D., Maignel, L., & Verrier, É. (1997). The value of using probabilities 
#' of gene origin to measure genetic variability in a population. Genetics Selection 
#' Evolution, 29(1), 5-23.
#'
#' @export
pedcontrib <- function(ped, reference = NULL, mode = c("both", "founder", "ancestor"), top = 20) {
  ped <- ensure_complete_tidyped(ped, "pedcontrib()")
  
  # Inject phantom parents to conserve genetic contributions accurately
  ped <- fill_phantoms(ped)
  
  mode <- match.arg(mode)
  
  # Ensure integer indices exist
  if (!all(c("SireNum", "DamNum", "IndNum") %in% names(ped))) {
    ped <- tidyped(ped, addnum = TRUE)
  }
  
  # Define reference cohort
  if (is.null(reference)) {
    max_gen <- max(ped$Gen)
    reference <- ped[Gen == max_gen, Ind]
    message(sprintf("Using %d individuals from generation %d as reference population.",
                    length(reference), max_gen))
  } else {
    if (!all(reference %in% ped$Ind)) {
      missing <- reference[!reference %in% ped$Ind]
      warning(sprintf("Reference IDs not found in pedigree: %s", paste(head(missing, 5), collapse = ", ")))
    }
    reference <- reference[reference %in% ped$Ind]
  }
  
  if (length(reference) == 0) {
    stop("No valid reference individuals specified.")
  }
  
  n <- nrow(ped)
  n_ref <- length(reference)
  
  # Pre-compute integer arrays
  sire_num <- ped$SireNum   # 0 = unknown
  dam_num <- ped$DamNum     # 0 = unknown
  ind_ids <- ped$Ind
  
  # Map reference IDs to integer positions
  ind_to_pos <- setNames(seq_len(n), ind_ids)
  cohort_pos <- as.integer(ind_to_pos[reference])
  
  # Determine numeric mode for C++: 1=founder, 2=ancestor, 3=both
  cpp_mode <- switch(mode, founder = 1L, ancestor = 2L, both = 3L)
  
  # ---- Call Rcpp backend ----
  cpp_res <- cpp_pedcontrib(
    sire     = as.integer(sire_num),
    dam      = as.integer(dam_num),
    cohort_pos = cohort_pos,
    mode     = cpp_mode
  )
  
  result <- list()
  
  # ---- Founder Contributions ----
  fe2_sum <- NA_real_
  fe_H <- NA_real_
  n_founders_total <- 0L
  
  if (mode %in% c("founder", "both")) {
    message("Calculating founder contributions...")
    
    f_idx <- cpp_res$founder_idx
    f_val <- cpp_res$founder_contrib
    
    if (length(f_idx) == 0) {
      warning("No founders found (all individuals have at least one parent).")
      result$founders <- data.table(Ind = character(0), Contrib = numeric(0),
                                     CumContrib = numeric(0), Rank = integer(0))
    } else {
      founder_dt <- data.table(
        Ind = ind_ids[f_idx],
        Contrib = f_val
      )
      
      setorder(founder_dt, -Contrib)
      founder_dt[, CumContrib := cumsum(Contrib)]
      founder_dt[, Rank := seq_len(.N)]
      
      fe2_sum <- sum(founder_dt$Contrib^2, na.rm = TRUE)
      fe_H <- calc_shannon_effective(founder_dt$Contrib)
      n_founders_total <- nrow(founder_dt)
      
      result$founders <- if (nrow(founder_dt) > top) founder_dt[1:top] else founder_dt
    }
  }
  
  # ---- Ancestor Contributions ----
  fa2_sum <- NA_real_
  fa_H <- NA_real_
  n_ancestors_total <- 0L
  
  if (mode %in% c("ancestor", "both")) {
    message("Calculating ancestor contributions (Boichard's iterative algorithm)...")
    
    a_idx <- cpp_res$ancestor_idx
    a_val <- cpp_res$ancestor_contrib
    
    if (length(a_idx) == 0) {
      result$ancestors <- data.table(Ind = character(0), Contrib = numeric(0),
                                      CumContrib = numeric(0), Rank = integer(0))
    } else {
      ancestor_dt_full <- data.table(
        Ind = ind_ids[a_idx],
        Contrib = a_val
      )
      ancestor_dt_full[, CumContrib := cumsum(Contrib)]
      ancestor_dt_full[, Rank := seq_len(.N)]
      
      fa2_sum <- sum(ancestor_dt_full$Contrib^2, na.rm = TRUE)
      fa_H <- calc_shannon_effective(ancestor_dt_full$Contrib)
      n_ancestors_total <- nrow(ancestor_dt_full)
      
      result$ancestors <- if (nrow(ancestor_dt_full) > top) ancestor_dt_full[1:top] else ancestor_dt_full
    }
  }
  
  # ---- Summary Statistics ----
  summary_list <- list(n_ref = n_ref)
  
  if (!is.null(result$founders)) {
    summary_list$n_founder <- n_founders_total
    summary_list$n_founder_show <- nrow(result$founders)
    summary_list$f_e <- if (!is.na(fe2_sum) && fe2_sum > 0) 1 / fe2_sum else NA_real_
    summary_list$f_e_H <- fe_H
  }
  
  if (!is.null(result$ancestors)) {
    summary_list$n_ancestor <- n_ancestors_total
    summary_list$n_ancestor_show <- nrow(result$ancestors)
    summary_list$f_a <- if (!is.na(fa2_sum) && fa2_sum > 0) 1 / fa2_sum else NA_real_
    summary_list$f_a_H <- fa_H
  }
  
  result$summary <- summary_list
  class(result) <- "pedcontrib"
  
  return(result)
}

#' Calculate Ancestry Proportions
#'
#' Estimates the proportion of genes for each individual that originates from 
#' specific founder groups (e.g., breeds, source populations).
#'
#' @param ped A \code{tidyped} object.
#' @param foundervar Character. The name of the column containing founder-group
#'   labels (e.g., "Breed", "Origin").
#' @param target_labels Character vector. Specific founder-group labels to track.
#'   If NULL, all unique labels in \code{foundervar} among founders are used.
#'
#' @return A \code{data.table} with columns:
#' \itemize{
#'   \item \code{Ind}: Individual ID.
#'   \item One column per tracked label (named after each unique value in
#'     \code{foundervar} among founders, or as specified by
#'     \code{target_labels}).
#'     Each value gives the proportion of genes (0--1) originating from that
#'     founder group. Row sums across all label columns equal 1.
#' }
#'
#' @examples
#' \donttest{
#' library(data.table)
#' # Create dummy labels for founders
#' tp <- tidyped(small_ped)
#' tp_dated <- copy(tp)
#' founders <- tp_dated[is.na(Sire) & is.na(Dam), Ind]
#' # Assign 'LineA' and 'LineB'
#' tp_dated[Ind %in% founders[1:(length(founders)/2)], Origin := "LineA"]
#' tp_dated[is.na(Origin), Origin := "LineB"]
#' 
#' # Calculate ancestry proportions for all individuals
#' anc <- pedancestry(tp_dated, foundervar = "Origin")
#' print(tail(anc))
#' }
#' 
#' @export
pedancestry <- function(ped, foundervar, target_labels = NULL) {
  ped <- ensure_complete_tidyped(ped, "pedancestry()")
  if (!foundervar %in% names(ped)) stop(sprintf("Column '%s' not found.", foundervar))
  
  # Forward pass requires numbers
  if (!all(c("SireNum", "DamNum", "IndNum") %in% names(ped))) {
    ped_work <- tidyped(ped, addnum = TRUE)
  } else {
    ped_work <- copy(ped)
  }
  
  # Ensure sorted by Gen for correct forward propagation
  setorder(ped_work, Gen, IndNum)
  
  # Identify founder labels
  founders_idx <- which(is.na(ped_work$Sire) & is.na(ped_work$Dam))
  if (is.null(target_labels)) {
    target_labels <- unique(ped_work[[foundervar]][founders_idx])
    target_labels <- target_labels[!is.na(target_labels) & target_labels != ""]
  }
  
  if (length(target_labels) == 0) {
    stop("No valid founder-group labels found. Please check 'foundervar'.")
  }

  n <- nrow(ped_work)
  res_mat <- matrix(0, nrow = n, ncol = length(target_labels))
  colnames(res_mat) <- target_labels
  
  # Initialize founders with their respective labels using vectorized matrix indexing
  founder_labels <- ped_work[[foundervar]][founders_idx]
  valid_matches <- match(founder_labels, target_labels)
  valid_mask <- !is.na(valid_matches)
  
  if (any(valid_mask)) {
    res_mat[cbind(founders_idx[valid_mask], valid_matches[valid_mask])] <- 1.0
  }
  
  sire_nums <- ped_work$SireNum
  dam_nums <- ped_work$DamNum
  
  # Forward pass via Rcpp
  res_mat <- cpp_calculate_ancestry(
    sire_nums, 
    dam_nums, 
    res_mat
  )
  
  ans <- as.data.table(res_mat)
  colnames(ans) <- target_labels
  ans[, Ind := ped_work$Ind]
  setcolorder(ans, "Ind")
  return(ans[])
}

#' Summarize Inbreeding Levels
#'
#' Classifies individuals into inbreeding levels based on their inbreeding 
#' coefficients (F) according to standard or user-defined thresholds.
#'
#' @param ped A \code{tidyped} object.
#' @param breaks Numeric vector of strictly increasing positive upper bounds for
#'   inbreeding classes. Default is \code{c(0.0625, 0.125, 0.25)}, corresponding
#'   approximately to half-sib, avuncular/grandparent, and full-sib/parent-offspring
#'   mating thresholds. The class \code{"F = 0"} is always kept as a fixed first
#'   level. A final open-ended class \code{"F > max(breaks)"} is always appended
#'   automatically.
#' @param labels Optional character vector of interval labels. If \code{NULL},
#'   labels are generated automatically from \code{breaks}. When supplied, its
#'   length must equal \code{length(breaks)}, with each element naming the
#'   bounded interval \code{(breaks[i-1], breaks[i]]}. The open-ended tail
#'   class is always auto-generated and cannot be overridden.
#'
#' @return A \code{data.table} with 3 columns:
#' \describe{
#'   \item{\code{FClass}}{An ordered factor. By default it contains 5 levels:
#'   \code{"F = 0"}, \code{"0 < F <= 0.0625"}, \code{"0.0625 < F <= 0.125"},
#'   \code{"0.125 < F <= 0.25"}, and \code{"F > 0.25"}. The number of levels
#'   equals \code{length(breaks) + 2} (the fixed zero class plus one class per
#'   bounded interval plus the open-ended tail).}
#'   \item{\code{Count}}{Integer. Number of individuals in each class.}
#'   \item{\code{Percentage}}{Numeric. Percentage of individuals in each class.}
#' }
#'
#' @details
#' The default thresholds follow common pedigree interpretation rules:
#' \itemize{
#'   \item \code{F = 0.0625}: approximately the offspring of half-sib mating.
#'   \item \code{F = 0.125}: approximately the offspring of avuncular or grandparent-grandchild mating.
#'   \item \code{F = 0.25}: approximately the offspring of full-sib or parent-offspring mating.
#' }
#' Therefore, assigning \code{F = 0.25} to the class \code{"0.125 < F <= 0.25"}
#' is appropriate. If finer reporting is needed, supply custom \code{breaks},
#' for example to separate \code{0.25}, \code{0.375}, or \code{0.5}.
#'
#' @examples
#' tp <- tidyped(simple_ped, addnum = TRUE)
#' pedfclass(tp)
#'
#' # Finer custom classes (4 breaks, labels auto-generated)
#' pedfclass(tp, breaks = c(0.03125, 0.0625, 0.125, 0.25))
#'
#' # Custom labels aligned to breaks (3 labels for 3 breaks; tail is auto)
#' pedfclass(tp, labels = c("Low", "Moderate", "High"))
#'
#' \donttest{
#' tp_inbred <- tidyped(inbred_ped, addnum = TRUE)
#' pedfclass(tp_inbred)
#' }
#'
#' @export
pedfclass <- function(ped,
                      breaks = c(0.0625, 0.125, 0.25),
                      labels = NULL) {
  ped <- ensure_tidyped(ped)

  if (!is.numeric(breaks) || length(breaks) < 1L || any(!is.finite(breaks)) ||
      any(breaks <= 0) || is.unsorted(breaks, strictly = TRUE)) {
    stop("breaks must be a strictly increasing numeric vector of positive values")
  }

  # Build bounded-interval labels (one per break, no tail)
  if (is.null(labels)) {
    lower_bounds <- c(0, head(breaks, -1))
    labels <- ifelse(
      lower_bounds == 0,
      paste0("0 < F <= ", breaks),
      paste0(lower_bounds, " < F <= ", breaks)
    )
  } else {
    if (!is.character(labels) || length(labels) != length(breaks)) {
      stop("labels must be a character vector of length equal to length(breaks)")
    }
  }
  # The open-ended tail class is always auto-generated
  tail_label <- paste0("F > ", max(breaks))
  all_labels  <- c(labels, tail_label)
  
  if (!"f" %in% names(ped)) {
    ped <- ensure_complete_tidyped(ped, "pedfclass()")
    message("Calculating inbreeding coefficients...")
    ped <- inbreed(ped)
  }
  
  cut_breaks   <- c(-Inf, .Machine$double.eps, breaks, Inf)
  class_levels <- c("F = 0", all_labels)

  # Classify
  f_classes <- cut(ped$f, breaks = cut_breaks, labels = class_levels, include.lowest = TRUE)
  
  dt <- data.table(FClass = factor(f_classes, levels = class_levels, ordered = TRUE))
  summary_dt <- dt[, .(Count = .N), by = FClass]
  
  # Ensure all categories appear
  all_cats <- data.table(FClass = factor(class_levels, levels = class_levels, ordered = TRUE))
  summary_dt <- merge(all_cats, summary_dt, by = "FClass", all.x = TRUE)
  summary_dt[is.na(Count), Count := 0]
  summary_dt[, Count := as.integer(Count)]
  
  summary_dt[, Percentage := (Count / sum(Count)) * 100]
  setorder(summary_dt, FClass)
  
  return(summary_dt[])
}

#' Calculate Partial Inbreeding
#'
#' Decomposes individuals' inbreeding coefficients into marginal contributions from 
#' specific ancestors. This allows identifying which ancestors or lineages are 
#' responsible for the observed inbreeding.
#'
#' @param ped A \code{tidyped} object.
#' @param ancestors Character vector. IDs of ancestors to calculate partial 
#'   inbreeding for. If NULL, the top ancestors by marginal contribution are used.
#' @param top Integer. Number of top ancestors to include if \code{ancestors} is NULL.
#'
#' @return A \code{data.table} with the first column as \code{Ind} and subsequent 
#'   columns representing the partial inbreeding ($pF$) from each ancestor.
#'
#' @details
#' The sum of all partial inbreeding coefficients for an individual (including 
#' contributions from founders) equals $1 + f_i$, where $f_i$ is the total 
#' inbreeding coefficient. This function specifically isolates the terms in 
#' the Meuwissen & Luo (1992) decomposition that correspond to the selected ancestors.
#'
#' @references
#' Lacey, R. C. (1996). A formula for determining the partial inbreeding coefficient, 
#' \eqn{F_{ij}}. Journal of Heredity, 87(4), 337-339.
#' 
#' Meuwissen, T. H., & Luo, Z. (1992). Computing inbreeding coefficients in 
#' large populations. Genetics Selection Evolution, 24(4), 305-313.
#'
#' @examples
#' \donttest{
#' library(data.table)
#' tp <- tidyped(inbred_ped)
#' # Calculate partial inbreeding originating from specific ancestors
#' target_ancestors <- inbred_ped[is.na(Sire) & is.na(Dam), Ind]
#' pF <- pedpartial(tp, ancestors = target_ancestors)
#' print(tail(pF))
#' }
#' 
#' @export
pedpartial <- function(ped, ancestors = NULL, top = 20) {
  ped <- ensure_complete_tidyped(ped, "pedpartial()")
  
  # Need dii from inbreeding algorithm which requires integer indices
  if (!all(c("SireNum", "DamNum", "IndNum") %in% names(ped))) {
    ped_work <- tidyped(ped, addnum = TRUE)
  } else {
    ped_work <- copy(ped)
  }
  
  # Ensure sorted by IndNum for alignment with C++ output
  setorder(ped_work, IndNum)
  
  # Identify target ancestors
  if (is.null(ancestors)) {
    message("No ancestors specified. Identifying top contributors...")
    cont <- pedcontrib(ped, mode = "ancestor", top = top)
    ancestors <- cont$ancestors$Ind
    if (length(ancestors) == 0) stop("Could not identify top ancestors.")
  }
  
  # Map ancestor IDs to IndNum
  anc_nums <- ped_work[Ind %in% ancestors, IndNum]
  anc_ids <- ped_work[Ind %in% ancestors, Ind]
  
  if (length(anc_nums) == 0) stop("None of the specified ancestors found in pedigree.")
  
  message(sprintf("Calculating partial inbreeding for %d ancestors...", length(anc_nums)))
  
  # Calculate dii (Meuwissen & Luo engines)
  # We reuse the logic from cpp_calculate_inbreeding but only need dii
  res_f <- cpp_calculate_inbreeding(ped_work$SireNum, ped_work$DamNum)
  
  # Adjust dii for half-founders to match Gulisija & Crow (2007) tabular method.
  # For half-founders (one parent known, one unknown), the standard dii includes

  # 0.25 Mendelian sampling variance from the unknown parent that cannot be
  # attributed to any known ancestor. We subtract this to ensure that
  # sum(pF_j) over all founder ancestors j equals the total F_i.
  # Standard dii for half-founder: 0.75 - 0.25 * F_known
  # Adjusted dii for partial F:    0.50 - 0.25 * F_known
  dii_partial <- res_f$dii
  sire_vec <- ped_work$SireNum
  dam_vec <- ped_work$DamNum
  f_vec <- res_f$f
  half_founders <- which(xor(sire_vec > 0, dam_vec > 0))
  if (length(half_founders) > 0) {
    for (idx in half_founders) {
      known_parent <- if (sire_vec[idx] > 0) sire_vec[idx] else dam_vec[idx]
      dii_partial[idx] <- 0.5 - 0.25 * f_vec[known_parent]
    }
  }
  
  # Call C++ partial inbreeding engine
  pf_mat <- cpp_calculate_partial_inbreeding(
    ped_work$SireNum, 
    ped_work$DamNum, 
    dii_partial, 
    as.integer(anc_nums)
  )
  
  result <- as.data.table(pf_mat)
  colnames(result) <- anc_ids
  result[, Ind := ped_work$Ind]
  setcolorder(result, "Ind")
  
  return(result[])
}

#' Print Founder and Ancestor Contributions
#'
#' @param x A \code{pedcontrib} object.
#' @param ... Additional arguments.
#'
#' @export
print.pedcontrib <- function(x, ...) {
  cat("Founder and Ancestor Contributions\n")
  cat("===================================\n")
  
  s <- x$summary
  cat(sprintf("Reference population size: %d\n", s$n_ref))
  
  if (!is.null(x$founders)) {
    cat(sprintf("\nFounders: %d (reported top %d)\n", s$n_founder, s$n_founder_show))
    cat(sprintf("  f_e(H) = %.3f  |  f_e = %.3f\n",
                if (is.null(s$f_e_H) || is.na(s$f_e_H)) NA_real_ else s$f_e_H,
                s$f_e))
    cat("\nTop 10 Founder Contributions:\n")
    print(head(x$founders, 10))
  }
  
  if (!is.null(x$ancestors)) {
    cat(sprintf("\nAncestors: %d (reported top %d)\n", s$n_ancestor, s$n_ancestor_show))
    cat(sprintf("  f_a(H) = %.3f  |  f_a = %.3f\n",
                if (is.null(s$f_a_H) || is.na(s$f_a_H)) NA_real_ else s$f_a_H,
                s$f_a))
    cat("\nTop 10 Ancestor Contributions:\n")
    print(head(x$ancestors, 10))
  }
  
  invisible(x)
}

#' Calculate Genetic Diversity Indicators
#'
#' Combines founder/ancestor contributions ($f_e$, $f_a$) and effective population
#' size estimates (Ne) from three methods into a single summary object.
#'
#' @param ped A \code{tidyped} object.
#' @param reference Character vector. Optional subset of individual IDs defining the reference population.
#'   If NULL, uses all individuals in the most recent generation.
#' @param top Integer. Number of top contributors to return in founder/ancestor tables. Default is 20.
#' @param nsamples Integer. Number of individuals sampled per cohort for the coancestry
#'   Ne method and for \eqn{f_g} estimation. Very large cohorts are sampled down to this
#'   size to control memory usage (default: 1000).
#' @param ncores Integer. Number of cores for parallel processing in the coancestry method. Default is 1.
#' @param seed Integer or NULL. Random seed passed to \code{set.seed()} before sampling
#'   in the coancestry method, ensuring reproducible \eqn{f_g} and \eqn{N_e} estimates.
#'   Default is \code{NULL} (no fixed seed).
#'
#' @return A list with class \code{pediv} containing:
#' \itemize{
#'   \item \code{summary}: A single-row \code{data.table} with columns
#'     \code{NRef}, \code{NFounder}, \code{feH}, \code{fe}, \code{NAncestor},
#'     \code{faH}, \code{fa}, \code{fafe}, \code{fg}, \code{MeanCoan},
#'     \code{NSampledCoan}, \code{NeCoancestry}, \code{NeInbreeding},
#'     \code{NeDemographic}.
#'     Here \code{feH} and \code{faH} are the Shannon-entropy (\eqn{q=1})
#'     effective numbers of founders and ancestors, respectively.
#'   \item \code{founders}: A \code{data.table} of top founder contributions.
#'   \item \code{ancestors}: A \code{data.table} of top ancestor contributions.
#' }
#'
#' @details
#' Internally calls \code{\link{pedcontrib}} for \eqn{f_e} and \eqn{f_a}.
#' The coancestry method is called via the internal \code{calc_ne_coancestry()}
#' function directly so that \eqn{f_g} and the Ne estimate can be obtained from
#' the same traced pedigree without duplication.
#' The inbreeding and demographic Ne methods are obtained via \code{\link{pedne}}.
#' All calculations use the same \code{reference} population.
#' If any method fails (e.g., insufficient pedigree depth), its value is \code{NA}
#' rather than stopping execution.
#'
#' \eqn{f_g} (founder genome equivalents, Caballero & Toro 2000) is estimated from
#' the diagonal-corrected mean coancestry of the reference population:
#' \deqn{\hat{\bar{C}} = \frac{N-1}{N} \cdot \frac{\bar{a}_{off}}{2} + \frac{1 + \bar{F}_s}{2N}}
#' \deqn{f_g = \frac{1}{2 \hat{\bar{C}}}}
#' where \eqn{N} is the full reference cohort size, \eqn{\bar{a}_{off}} is the
#' off-diagonal mean relationship among sampled individuals, and \eqn{\bar{F}_s}
#' is their mean inbreeding coefficient.
#'
#' @examples
#' \donttest{
#' tp <- tidyped(small_ped)
#' div <- pediv(tp, reference = c("Z1", "Z2", "X", "Y"), seed = 42L)
#' print(div)
#'
#' # Access Shannon effective numbers from summary
#' div$summary$feH   # Shannon effective founders (q=1)
#' div$summary$faH   # Shannon effective ancestors (q=1)
#'
#' # Founder diversity profile: NFounder >= feH >= fe
#' with(div$summary, c(NFounder = NFounder, feH = feH, fe = fe))
#' }
#'
#' @seealso \code{\link{pedcontrib}}, \code{\link{pedne}}, \code{\link{pedstats}}
#' @export
pediv <- function(ped, reference = NULL, top = 20, nsamples = 1000, ncores = 1,
                  seed = NULL) {
  ped <- ensure_complete_tidyped(ped, "pediv()")

  # ---- Founder and ancestor contributions ----
  message("Calculating founder and ancestor contributions...")
  contrib <- pedcontrib(ped, reference = reference, mode = "both", top = top)
  s <- contrib$summary

  # ---- Prepare pedigree with ECG (mirrors pedne() preprocessing) ----
  ped_dt <- data.table::copy(ped)
  if (!"ECG" %in% names(ped_dt)) {
    ecg_dt <- pedecg(ped_dt)
    ped_dt <- merge(ped_dt, ecg_dt[, .(Ind, ECG)], by = "Ind", all.x = TRUE)
  }
  by_all <- ".CohortAll"
  ped_dt[[by_all]] <- "All"

  if (!is.null(reference)) {
    ped_subset <- ped_dt[Ind %in% reference]
  } else {
    ped_subset <- ped_dt
  }

  # ---- Ne (coancestry) + fg: call internal function directly ----
  # pedne(method="coancestry") strips fg from its public output; bypass it here
  # so both Ne and fg are computed from the same single traced pedigree.
  message("Calculating Ne (coancestry) and fg...")
  raw_coan <- tryCatch(
    calc_ne_coancestry(ped_subset, ped_dt, by_all, nsamples, seed = seed),
    error = function(e) NULL
  )

  # ---- Ne: inbreeding and demographic (via public pedne API) ----
  message("Calculating Ne (inbreeding)...")
  ne_i <- tryCatch(
    suppressMessages(pedne(ped, reference = reference, by = NULL,
                           method = "inbreeding")),
    error = function(e) NULL
  )

  message("Calculating Ne (demographic)...")
  ne_d <- tryCatch(
    suppressMessages(pedne(ped, reference = reference, by = NULL,
                           method = "demographic")),
    error = function(e) NULL
  )

  extract_ne <- function(dt) {
    if (is.null(dt) || nrow(dt) == 0) NA_real_ else dt$Ne[1]
  }

  # Extract fg and MeanCoan from raw coancestry result
  fg_val      <- if (!is.null(raw_coan) && nrow(raw_coan) > 0 && "fg" %in% names(raw_coan))
                   raw_coan$fg[1] else NA_real_
  mean_coan   <- if (!is.null(raw_coan) && nrow(raw_coan) > 0 && "MeanCoan" %in% names(raw_coan))
                   raw_coan$MeanCoan[1] else NA_real_
  n_samp_coan <- if (!is.null(raw_coan) && nrow(raw_coan) > 0 && "NSampledCoan" %in% names(raw_coan))
                   as.integer(raw_coan$NSampledCoan[1]) else NA_integer_

  # ---- Assemble summary ----
  summary_dt <- data.table::data.table(
    NRef          = s$n_ref,
    NFounder      = s$n_founder,
    feH           = if (!is.null(s$f_e_H)) s$f_e_H else NA_real_,
    fe            = s$f_e,
    NAncestor     = s$n_ancestor,
    faH           = if (!is.null(s$f_a_H)) s$f_a_H else NA_real_,
    fa            = s$f_a,
    fafe          = if (!is.na(s$f_a) && !is.na(s$f_e) && s$f_e > 0)
                      round(s$f_a / s$f_e, 6) else NA_real_,
    fg            = fg_val,
    MeanCoan      = mean_coan,
    NSampledCoan  = n_samp_coan,
    NeCoancestry  = extract_ne(raw_coan),
    NeInbreeding  = extract_ne(ne_i),
    NeDemographic = extract_ne(ne_d)
  )

  result <- list(
    summary   = summary_dt,
    founders  = contrib$founders,
    ancestors = contrib$ancestors
  )
  class(result) <- "pediv"
  return(result)
}

#' Print Genetic Diversity Summary
#'
#' @param x A \code{pediv} object.
#' @param ... Additional arguments.
#'
#' @export
print.pediv <- function(x, ...) {
  cat("Genetic Diversity Summary\n")
  cat("=========================\n")

  s <- x$summary

  cat(sprintf("Reference population size : %d\n", s$NRef))

  cat("\n-- Founder / Ancestor Contributions --\n")
  cat(sprintf("Founders  (total) : %d    fe(H) = %.3f    fe = %.3f\n",
              s$NFounder, s$feH, s$fe))
  cat(sprintf("Ancestors (total) : %d    fa(H) = %.3f    fa = %.3f\n",
              s$NAncestor, s$faH, s$fa))
  cat(sprintf("fa/fe ratio       : %.4f\n", s$fafe))
  if (!is.null(s$fg) && !is.na(s$fg)) {
    cat(sprintf("Founder genomes   : fg = %.3f  (MeanCoan = %.6f, NSampled = %d)\n",
                s$fg, s$MeanCoan, s$NSampledCoan))
    # Hierarchy fg <= fa <= fe <= NFounder holds only when all metrics share
    # the same reference population. Display it only when the order is intact.
    if (!is.na(s$fa) && !is.na(s$fe) && s$fg <= s$fa + 1e-6 && s$fa <= s$fe + 1e-6) {
      cat(sprintf("Hierarchy: fg <= fa <= fe <= NFounder  =  %.2f <= %.3f <= %.3f <= %d\n",
                  s$fg, s$fa, s$fe, s$NFounder))
    }
  }

  cat("\n-- Effective Population Size (Ne) --\n")
  fmt_ne <- function(v) if (is.na(v)) "     NA" else sprintf("%7.1f", v)
  cat(sprintf("  Coancestry  : %s\n", fmt_ne(s$NeCoancestry)))
  cat(sprintf("  Inbreeding  : %s\n", fmt_ne(s$NeInbreeding)))
  cat(sprintf("  Demographic : %s\n", fmt_ne(s$NeDemographic)))

  if (!is.null(x$founders) && nrow(x$founders) > 0) {
    cat("\nTop Founder Contributions (top 5):\n")
    print(head(x$founders, 5))
  }

  if (!is.null(x$ancestors) && nrow(x$ancestors) > 0) {
    cat("\nTop Ancestor Contributions (top 5):\n")
    print(head(x$ancestors, 5))
  }

  invisible(x)
}
