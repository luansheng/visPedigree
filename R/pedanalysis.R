#' Calculate Generation Intervals
#'
#' Computes the generation intervals for the four gametic pathways: 
#' Sire to Son (SS), Sire to Daughter (SD), Dam to Son (DS), and Dam to Daughter (DD).
#' The generation interval is defined as the age of the parents at the birth of their offspring.
#'
#' @param ped A \code{tidyped} object.
#' @param timevar Character. The name of the column containing time information (e.g., "BirthYear", "Year").
#'   Defaults to "Year" or "BirthYear" if found.
#' @param unit Character. Time unit for the interval: \code{"year"} (default), 
#'   \code{"month"}, \code{"day"}, \code{"hour"}, or \code{"gen"} (raw numeric).
#' @param format Character. Optional format string for parsing \code{timevar} if it's character.
#' @param cycle_length Numeric. Optional target (designed) length of one generation
#'   cycle expressed in \code{unit}s. When provided, an additional column
#'   \code{GenEquiv} is appended to the result, defined as:
#'   \deqn{GenEquiv_i = \frac{\bar{L}_i}{L_{cycle}}}
#'   where \eqn{\bar{L}_i} is the observed mean interval for pathway \eqn{i} and
#'   \eqn{L_{cycle}} is \code{cycle_length}. A value > 1 means the observed
#'   interval exceeds the target cycle (lower breeding efficiency).
#'   Example: for Pacific white shrimp with a 180-day target cycle, set
#'   \code{unit = "day", cycle_length = 180}.
#' @param by Character. Optional grouping column (e.g., "Breed", "Farm"). 
#'   If provided, intervals are calculated within each group.
#'
#' @return A \code{data.table} with columns:
#' \itemize{
#'   \item \code{Group}: Grouping level (if \code{by} is used).
#'   \item \code{Pathway}: One of "SS", "SD", "DS", "DD", "SO", "DO", or "Average".
#'     SS/SD/DS/DD require offspring sex; SO (Sire-Offspring) and DO (Dam-Offspring)
#'     are computed from all parent-offspring pairs regardless of offspring sex.
#'   \item \code{N}: Number of parent-offspring pairs used.
#'   \item \code{Mean}: Average generation interval in \code{unit}.
#'   \item \code{SD}: Standard deviation of the interval.
#'   \item \code{GenEquiv}: (Optional) Generation equivalents based on \code{cycle_length}.
#' }
#' 
#' @examples
#' \dontrun{
#' # ---- Case 1: Integer or numeric year column (most common) ----
#' # Pedigree column 'BirthYear' contains values like 2020, 2021, 2022
#' tped <- tidyped(ped)   # ped must have a BirthYear column
#' pedgenint(tped, timevar = "BirthYear")  # unit defaults to "year"
#'
#' # ---- Case 2: Standard ISO date strings "YYYY-MM-DD" ----
#' # Pedigree column 'HatchDate' contains values like "2020-06-15"
#' # format= is NOT needed; the function detects ISO strings automatically.
#' pedgenint(tped, timevar = "HatchDate", unit = "day")
#'
#' # ---- Case 3: Non-standard date format, must specify format= ----
#' # Pedigree column 'HatchDate' contains values like "15/06/2020" (DD/MM/YYYY)
#' pedgenint(tped, timevar = "HatchDate", unit = "day", format = "%d/%m/%Y")
#'
#' # ---- Case 4: cycle_length - generation equivalents ----
#' # Pacific white shrimp: target one generation per 180 days.
#' # GenEquiv = Mean / 180; value > 1 means observed interval exceeds target.
#' pedgenint(tped, timevar = "HatchDate", unit = "day", cycle_length = 180)
#'
#' # ---- Case 5: Grouping by breed or farm ----
#' # Pedigree has an additional column 'Breed'
#' pedgenint(tped, timevar = "BirthYear", by = "Breed")
#' }
#'
#' @export
pedgenint <- function(ped, timevar = NULL, unit = c("year", "month", "day", "hour", "gen"), 
                      format = NULL, cycle_length = NULL, by = NULL) {
  # Input validation
  if (!inherits(ped, "tidyped")) stop("ped must be a tidyped object")
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
  
  # Convert ID to character for matching
  dt[, `:=`(Ind = as.character(Ind), 
            Sire = as.character(Sire), 
            Dam = as.character(Dam))]
  
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
  if (!is.null(cycle_length) && is.numeric(cycle_length)) {
    final_res[, GenEquiv := Mean / cycle_length]
  }
  
  if (!is.null(by)) {
    setorderv(final_res, c("Group", "Pathway"))
  } else {
    setorder(final_res, Pathway)
  }
  
  attr(final_res, "unit") <- unit
  return(final_res)
}

.parse_to_numeric_time <- function(x, unit = "year", format = NULL) {
  # Internal factors relative to second
  factors <- c(year = 31557600, month = 2629800, day = 86400, hour = 3600, gen = 1)
  target_factor <- factors[unit]
  
  if (is.numeric(x)) {
    # If numeric, we assume input is in 'years' by default unless unit is 'gen'
    if (unit == "gen") return(x)
    if (unit == "year") return(x)
    # Convert 'years' to target unit
    return(x * (factors["year"] / target_factor))
  }
  
  # String parsing with caching for performance
  u_x <- unique(x[!is.na(x)])
  if (length(u_x) == 0) return(rep(NA_real_, length(x)))
  
  parsed <- if (!is.null(format)) {
    as.POSIXct(u_x, format = format)
  } else {
    # Try ISO/Common formats
    p <- suppressWarnings(as.POSIXct(u_x, optional = TRUE))
    if (all(is.na(p))) {
      p <- suppressWarnings(as.Date(u_x, optional = TRUE))
    }
    if (all(is.na(p))) {
      # Fallback for YYYY/MM/DD
      p <- suppressWarnings(as.Date(u_x, format = "%Y/%m/%d", optional = TRUE))
    }
    as.POSIXct(p)
  }
  
  if (all(is.na(parsed))) {
    stop("Failed to parse time variable. Please provide numeric years or specify 'format' (e.g. '%Y-%m-%d').")
  }
  
  val_map <- setNames(as.numeric(parsed), u_x)
  numeric_seconds <- val_map[as.character(x)]
  
  return(numeric_seconds / target_factor)
}

#' Pedigree Subpopulations
#'
#' Splits the pedigree into subpopulations and calculates basic statistics for each group.
#'
#' @param ped A \code{tidyped} object.
#' @param by Character. The name of the column to group by. 
#'   If NULL, uses connected components via \code{\link{splitped}}.
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
#' @export
pedsubpop <- function(ped, by = NULL) {
  if (!inherits(ped, "tidyped")) stop("ped must be a tidyped object")
  
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
  
  data.table::rbindlist(res_list)
}

#' Calculate Average Relationship Coefficients
#'
#' Computes the average pairwise relationship coefficients within cohorts or groups.
#' This corresponds to the "mean relationship" mentioned in academic literature 
#' to monitor genetic diversity.
#'
#' @param ped A \code{tidyped} object.
#' @param by Character. The column name to group by (e.g., "Year", "Breed", "Generation").
#' @param cand Character vector. An optional vector of candidate individual IDs to calculate 
#'   relationships for. If provided, only individuals matching these IDs in each group 
#'   will be used. Default is NULL (use all individuals in the group).
#' @param compact Logical. Whether to use compact representation for large families to 
#'   save memory. Recommended when pedigree size exceeds 25,000. Default is FALSE.
#'
#' @return A \code{data.table} with columns:
#' \itemize{
#'   \item \code{Group}: The grouping identifier.
#'   \item \code{NTotal}: Total number of individuals in the group.
#'   \item \code{NUsed}: Number of individuals used in calculation (could be subset by cand).
#'   \item \code{MeanRel}: Average of off-diagonal elements in the A matrix for this group. NA if less than 2 individuals.
#' }
#' 
#' @export
pedrel <- function(ped, by = "Gen", cand = NULL, compact = FALSE) {
  if (!inherits(ped, "tidyped")) stop("ped must be a tidyped object")
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
    
    if (!is.null(cand)) {
      sub_ped <- sub_ped_full[Ind %in% cand]
    } else {
      sub_ped <- sub_ped_full
    }
    n_used <- nrow(sub_ped)
    
    if (n_used < 2) {
      warning(sprintf("Group '%s' has less than 2 individuals after applying 'cand', returning NA_real_.", g))
      return(data.table(Group = g, NTotal = n_total, NUsed = n_used, MeanRel = NA_real_))
    }
    
    local_ped <- tryCatch({
      suppressMessages(tidyped(sub_ped, addnum = TRUE))
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
                 "Hint: use 'cand' parameter to reduce group size, or set compact = TRUE."),
          g, nrow(local_ped), max_dense))
        return(data.table(Group = g, NTotal = n_total, NUsed = n_used, MeanRel = NA_real_))
      }
      
      if (isTRUE(compact)) {
        A <- tryCatch({
          pedmat(local_ped, method = "A", sparse = FALSE, compact = TRUE)
        }, error = function(e) {
          warning(sprintf(
            paste0("Group '%s': %s\n",
                   "Hint: subset with 'cand' or use 'compact = TRUE'."),
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
                   "Hint: subset with 'cand' or use 'compact = TRUE'."),
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
  return(result)
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
#' @references
#' Boichard, D., Maignel, L., & Verrier, E. (1997). The value of using probabilities of gene origin to measure genetic variability in a population. Genetics Selection Evolution, 29(1), 5.
#'
#' @export
pedecg <- function(ped) {
  if (!inherits(ped, "tidyped")) stop("ped must be a tidyped object")
  
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
    
    valid_s <- sires > 0
    if (any(valid_s)) {
      parents <- sires[valid_s]
      sire_ecg[valid_s] <- ecg[parents]
      sire_full[valid_s] <- full_gen[parents]
      sire_max[valid_s] <- max_gen_val[parents]
    }
    
    valid_d <- dams > 0
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
  
  data.table(
    Ind = ped$Ind,
    ECG = ecg,
    FullGen = full_gen,
    MaxGen = max_gen_val
  )
}

#' Pedigree Statistics
#'
#' Calculates comprehensive statistics for a pedigree, including population structure,
#' generation intervals, and ancestral depth.
#'
#' @param ped A \code{tidyped} object.
#' @param timevar Optional character. Name of the column containing birth time or
#'   year information (e.g., \code{"BirthYear"}, \code{"HatchDate"}).
#'   Accepted column formats:
#'   \itemize{
#'     \item \strong{Integer/numeric year}: e.g., \code{2020}, \code{2020.5}.
#'       Use \code{unit = "year"} (default).
#'     \item \strong{ISO date string}: e.g., \code{"2020-06-15"}.
#'       Detected automatically; set \code{unit} to \code{"day"} or \code{"month"}.
#'     \item \strong{Custom date format}: e.g., \code{"15/06/2020"} (DD/MM/YYYY).
#'       Must specify \code{format = "\%d/\%m/\%Y"} and pass via \code{...}.
#'   }
#'   If \code{NULL}, attempts auto-detection from common column names
#'   (\code{"BirthYear"}, \code{"Year"}, \code{"BirthDate"}, etc.).
#' @param unit Character. Time unit for reporting generation intervals:
#'   \code{"year"} (default), \code{"month"}, \code{"day"}, \code{"hour"},
#'   or \code{"gen"} (raw numeric). Must match the scale of \code{timevar};
#'   e.g., use \code{unit = "day"} when \code{timevar} holds date strings.
#' @param cycle_length Numeric. Optional target generation cycle length in
#'   \code{unit}s. When provided, \code{gen_intervals} will include a
#'   \code{GenEquiv} column (observed Mean / cycle_length). See
#'   \code{\link{pedgenint}} for details.
#' @param ... Additional arguments passed to \code{\link{pedgenint}},
#'   e.g., \code{format} for custom date parsing or \code{by} for grouping.
#'
#' @return An object of class \code{pedstats}, which is a list containing:
#' \itemize{
#'   \item \code{summary}: Basic summary statistics.
#'   \item \code{gen_intervals}: Generation intervals (\code{NULL} if no
#'     \code{timevar} is detected).
#'   \item \code{ecg}: Equi-Generate Coefficients and ancestral depth.
#' }
#'
#' @examples
#' \dontrun{
#' # ---- Standard annual pedigree ----
#' # 'BirthYear' column contains integer years like 2020, 2021
#' tped <- tidyped(ped)
#' pedstats(tped, timevar = "BirthYear")
#'
#' # ---- Short-cycle species (e.g. Pacific white shrimp) ----
#' # 'HatchDate' column contains ISO strings like "2020-06-15"
#' # target generation cycle = 180 days
#' pedstats(tped, timevar = "HatchDate", unit = "day", cycle_length = 180)
#'
#' # ---- Custom date format passed via ... ----
#' # 'HatchDate' contains "15/06/2020" (DD/MM/YYYY)
#' pedstats(tped, timevar = "HatchDate", unit = "day", format = "%d/%m/%Y")
#' }
#'
#' @export
pedstats <- function(ped, timevar = NULL, unit = "year", cycle_length = NULL, calc_ecg = TRUE, calc_genint = TRUE, ...) {
  if (!inherits(ped, "tidyped")) stop("ped must be a tidyped object")
  
  summ <- data.table(
    N = nrow(ped),
    N_Sire = data.table::uniqueN(ped$Sire, na.rm = TRUE),
    N_Dam = data.table::uniqueN(ped$Dam, na.rm = TRUE),
    N_Founder = sum(is.na(ped$Sire) & is.na(ped$Dam)),
    Max_Gen = max(ped$Gen)
  )
  
  ecg_dt <- NULL
  if (calc_ecg) {
    ecg_dt <- pedecg(ped)
  }
  
  gen_int <- NULL
  
  target_timevar <- timevar
  if (is.null(target_timevar)) {
    candidates <- c("BirthYear", "Year", "birth_year", "year", "BirthDate", "Date")
    match_col <- intersect(names(ped), candidates)
    if (length(match_col) > 0) target_timevar <- match_col[1]
  }
  
  if (calc_genint && !is.null(target_timevar) && target_timevar %in% names(ped)) {
    gen_int <- tryCatch({
      pedgenint(ped, timevar = target_timevar, unit = unit, cycle_length = cycle_length, ...)
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
  cat(sprintf("Sires: %d | Dams: %d | Founders: %d\n", s$N_Sire, s$N_Dam, s$N_Founder))
  cat(sprintf("Max Generation: %d\n", s$Max_Gen))
  
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
#' Calculates the effective population size (Ne) based on individual rate of inbreeding (ΔF).
#' Uses the method described by Gutiérrez et al. (2008, 2009).
#'
#' @param ped A \code{tidyped} object.
#' @param by Character. The name of the column used to define cohorts (e.g., "Year", "BirthYear").
#'   If NULL, attempts to auto-detect from common names.
#' @param cand Character vector. Optional subset of individual IDs defining the reference cohort.
#'   If NULL, uses all individuals.
#'
#' @return A \code{data.table} with columns:
#' \itemize{
#'   \item \code{Cohort}: Year or period identifier.
#'   \item \code{N}: Number of individuals in cohort.
#'   \item \code{MeanF}: Mean inbreeding coefficient.
#'   \item \code{DeltaF}: Individual increase in inbreeding per generation.
#'   \item \code{Ne}: Effective population size.
#' }
#' 
#' @details
#' The effective population size is calculated using one of three methods:
#' 
#' \itemize{
#'   \item \strong{"coancestry"} (Default): Based on the rate of coancestry between pairs of individuals
#'   (Cervantes et al. 2011).
#'   \deqn{\Delta c_{ij} = 1 - (1 - c_{ij})^{1/(\frac{ECG_i + ECG_j}{2})}}
#'   \deqn{N_e = \frac{1}{2 \overline{\Delta c}}}
#'   To handle large populations, this method samples \code{nsamples} individuals per cohort
#'   and computes the mean rate of coancestry among them.
#'   
#'   \item \strong{"inbreeding"}: Based on the individual rate of inbreeding 
#'   (Gutiérrez et al. 2008, 2009).
#'   \deqn{\Delta F_i = 1 - (1 - F_i)^{1/(ECG_i - 1)}}
#'   \deqn{N_e = \frac{1}{2 \overline{\Delta F}}}
#'   
#'   \item \strong{"demographic"}: Based on the number of breeding males and females.
#'   \deqn{N_e = \frac{4 N_m N_f}{N_m + N_f}}
#'   Where \eqn{N_m} and \eqn{N_f} are the number of unique male and female parents
#'   contributing to the cohort.
#' }
#'
#' @param ped A tidyped object.
#' @param method Character, method to use: "coancestry" (default), "inbreeding", or "demographic".
#' @param by Character, column name for grouping cohorts (e.g., "Year"). If NULL, calculates overall Ne for all individuals.
#' @param cand Character vector, specific candidates to include.
#' @param nsamples Integer. Number of individuals to sample per cohort. Only applicable for \code{method = "coancestry"} (default: 1000).
#' @param ncores Integer. Number of cores for parallel processing. Currently only effective for \code{method = "coancestry"} (default: 1).
#'
#' @return A data.table with columns:
#' \itemize{
#'   \item \code{Cohort}: Grouping variable value
#'   \item \code{N}: Number of individuals
#'   \item \code{Ne}: Effective population size
#'   \item Additional columns depending on method (e.g., MeanF, DeltaF, Nm, Nf)
#' }
#' 
#' @export
pedne <- function(ped, method = c("coancestry", "inbreeding", "demographic"), 
                  by = NULL, cand = NULL, nsamples = 1000, ncores = 1) {
  
  if (!inherits(ped, "tidyped")) stop("ped must be a tidyped object")
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
  
  # Filter candidates if specified
  if (!is.null(cand)) {
    if (!all(cand %in% ped_dt$Ind)) {
      missing <- cand[!cand %in% ped_dt$Ind]
      warning(sprintf("Candidate IDs not found in pedigree: %s", paste(missing, collapse = ", ")))
    }
    ped_subset <- ped_dt[Ind %in% cand]
  } else {
    ped_subset <- ped_dt
  }
  
  # Dispatch to specific calculator
  res <- switch(method,
    "inbreeding" = calc_ne_inbreeding(ped_subset, by),
    "coancestry" = calc_ne_coancestry(ped_subset, by, nsamples),
    "demographic" = calc_ne_demographic(ped_subset, by)
  )
  
  return(res)
}

calc_ne_inbreeding <- function(ped, by) {
  # Ensure inbreeding coefficients and ECG are present
  if (!"f" %in% names(ped)) ped <- inbreed(ped)
  if (!"ECG" %in% names(ped)) {
    ecg_dt <- pedecg(ped)
    ped <- merge(ped, ecg_dt[, .(Ind, ECG)], by = "Ind", all.x = TRUE)
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
  
  return(result)
}

calc_ne_demographic <- function(ped, by) {
  # Ensure numbers are present
  if (!all(c("Sire", "Dam", "Sex") %in% names(ped))) {
     stop("Pedigree must contain 'Sire', 'Dam', and 'Sex' columns for demographic Ne.")
  }

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
  return(result)
}

calc_ne_coancestry <- function(ped, by, nsamples) {
  # Needs ECG for Delta C calculation
  if (!"ECG" %in% names(ped)) {
    ecg_dt <- pedecg(ped)
    ped <- merge(ped, ecg_dt[, .(Ind, ECG)], by = "Ind", all.x = TRUE)
  }
  
  ped[, CohortLabel := get(by)]
  cohorts <- sort(unique(ped$CohortLabel))
  
  results_list <- list()
  
  for (coh in cohorts) {
    sub_ped <- ped[CohortLabel == coh]
    n_total <- nrow(sub_ped)
    
    if (n_total == 0) next
    
    # Sampling
    if (n_total > nsamples) {
      sampled_inds <- sample(sub_ped$Ind, nsamples)
    } else {
      sampled_inds <- sub_ped$Ind
    }
    
    # Trace back ancestors to build sub-pedigree
    # This is critical for performance to avoid large matrix on full pedigree
    trace_ped <- tidyped(ped, cand = sampled_inds, trace = "up", addgen = TRUE, addnum = TRUE)
    
    # Ensure sorted by Gen for correct matrix build (though C++ uses numeric)
    setorder(trace_ped, Gen, IndNum)
    
    # Re-identify target numeric IDs
    target_nums <- trace_ped[Ind %in% sampled_inds, IndNum]
    
    # Ensure ECG is present
    if (!"ECG" %in% names(trace_ped)) {
       ecg_map <- ped[Ind %in% trace_ped$Ind, .(Ind, ECG)]
       trace_ped <- merge(trace_ped, ecg_map, by = "Ind", all.x = TRUE)
    }
    
    # If any ECG missing (e.g. base founders in trace not in original?), fill 0
    trace_ped[is.na(ECG), ECG := 0]
    
    # Call C++
    # Ensure vectors are 0-indexed for C++ if needed, but Rcpp usually handles R 1-based logic 
    # if implemented. Our C++ uses 1-based input (subtracts 1 inside).
    
    # The trace_ped MUST be reordered by IndNum 1..N
    setorder(trace_ped, IndNum)
    
    avg_delta_c <- tryCatch({
       cpp_calculate_sampled_coancestry_delta(
        trace_ped$SireNum, 
        trace_ped$DamNum, 
        target_nums, 
        trace_ped$ECG
      )
    }, error = function(e) {
      warning(paste("Error in coancestry calculation for cohort", coh, ":", e$message))
      return(NA_real_)
    })

    ne_val <- if (!is.na(avg_delta_c) && avg_delta_c > 0) 1 / (2 * avg_delta_c) else NA_real_
    
    results_list[[as.character(coh)]] <- data.table(
      Cohort = coh,
      N = n_total,
      NSampled = length(sampled_inds),
      DeltaC = avg_delta_c,
      Ne = ne_val
    )
  }
  
  result <- rbindlist(results_list)
  return(result)
}

# Legacy wrapper/stub for documentation or compatibility if needed
# (none required as we replaced pedne in place) -> THIS LINE INTENTIONALLY LEFT BLANK
#' Calculate Founder and Ancestor Contributions
#'
#' Calculates genetic contributions from founders (fe) and influential ancestors (fa).
#' Implements the gene dropping algorithm for founder contributions and Boichard's 
#' algorithm for ancestor contributions.
#'
#' @param ped A \code{tidyped} object.
#' @param cohort Character vector. Optional subset of individual IDs defining the reference population.
#'   If NULL, uses all individuals in the most recent generation.
#' @param mode Character. Type of contribution to calculate:
#' \itemize{
#'   \item \code{"founder"}: Founder contributions (fe).
#'   \item \code{"ancestor"}: Ancestor contributions (fa).
#'   \item \code{"both"}: Both founder and ancestor contributions.
#' }
#' @param top Integer. Number of top contributors to return. Default is 20.
#'
#' @return A list with class \code{pedcontrib} containing:
#' \itemize{
#'   \item \code{founders}: \code{data.table} of founder contributions (if mode includes "founder").
#'   \item \code{ancestors}: \code{data.table} of ancestor contributions (if mode includes "ancestor").
#'   \item \code{summary}: Summary statistics including effective number of founders/ancestors.
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
#' **Founder Contributions (fe):**
#' Calculated by probabilistic gene flow from founders to the reference cohort.
#' A founder is an individual with unknown parents.
#' 
#' **Ancestor Contributions (fa):**
#' Calculated using Boichard's iterative algorithm, accounting for:
#' \itemize{
#'   \item Marginal genetic contribution of each ancestor
#'   \item Long-term contributions through multiple pathways
#' }
#' 
#' **Effective Numbers:**
#' \itemize{
#'   \item \code{Ne_f}: Effective number of founders = 1 / sum(fe^2)
#'   \item \code{Ne_a}: Effective number of ancestors = 1 / sum(fa^2)
#' }
#' 
#' @references
#' Boichard, D., Maignel, L., & Verrier, É. (1997). The value of using probabilities 
#' of gene origin to measure genetic variability in a population. Genetics Selection 
#' Evolution, 29(1), 5-23.
#'
#' @export
pedcontrib <- function(ped, cand = NULL, mode = c("both", "founder", "ancestor"), top = 20) {
  if (!inherits(ped, "tidyped")) stop("ped must be a tidyped object")
  
  mode <- match.arg(mode)
  
  # Ensure integer indices exist
  if (!all(c("SireNum", "DamNum", "IndNum") %in% names(ped))) {
    ped <- tidyped(ped, addnum = TRUE)
  }
  
  # Define reference cohort
  if (is.null(cand)) {
    max_gen <- max(ped$Gen)
    cand <- ped[Gen == max_gen, Ind]
    message(sprintf("Using %d individuals from generation %d as reference candidate.",
                    length(cand), max_gen))
  } else {
    if (!all(cand %in% ped$Ind)) {
      missing <- cand[!cand %in% ped$Ind]
      warning(sprintf("Candidate IDs not found in pedigree: %s", paste(head(missing, 5), collapse = ", ")))
    }
    cand <- cand[cand %in% ped$Ind]
  }
  
  if (length(cand) == 0) {
    stop("No valid candidate individuals specified.")
  }
  
  n <- nrow(ped)
  n_cohort <- length(cand)
  
  # Pre-compute integer arrays for speed (avoid character matching)
  sire_num <- ped$SireNum   # 0 = unknown
  dam_num <- ped$DamNum     # 0 = unknown
  ind_ids <- ped$Ind
  
  # Map candidate IDs to integer positions
  ind_to_pos <- setNames(seq_len(n), ind_ids)
  cohort_pos <- ind_to_pos[cand]
  # Vectorized backward pass: compute total flow from all cohort members at once
  # contrib[i] = average probability that a randomly chosen gene from the cohort
  #              traces back to individual i
  backward_pass <- function(sire_v, dam_v, cohort_positions, n_ped, n_coh) {
    contrib <- numeric(n_ped)
    contrib[cohort_positions] <- 1.0 / n_coh
    
    # Walk backwards (individuals sorted by Gen ascending, so reverse = top-down)
    for (i in n_ped:1) {
      if (contrib[i] == 0) next
      s <- sire_v[i]
      d <- dam_v[i]
      if (s > 0) contrib[s] <- contrib[s] + 0.5 * contrib[i]
      if (d > 0) contrib[d] <- contrib[d] + 0.5 * contrib[i]
    }
    return(contrib)
  }
  
  result <- list()
  
  # ---- Founder Contributions ----
  fe2_sum <- NA_real_
  n_founders_total <- 0L
  
  if (mode %in% c("founder", "both")) {
    message("Calculating founder contributions...")
    
    founder_idx <- which(sire_num == 0 & dam_num == 0)
    n_founders <- length(founder_idx)
    
    if (n_founders == 0) {
      warning("No founders found (all individuals have at least one parent).")
      result$founders <- data.table(Ind = character(0), Contrib = numeric(0),
                                     CumContrib = numeric(0), Rank = integer(0))
    } else {
      # Single vectorized backward pass for all cohort members simultaneously
      contrib <- backward_pass(sire_num, dam_num, cohort_pos, n, n_cohort)
      
      founder_dt <- data.table(
        Ind = ind_ids[founder_idx],
        Contrib = contrib[founder_idx]
      )
      
      setorder(founder_dt, -Contrib)
      founder_dt[, CumContrib := cumsum(Contrib)]
      founder_dt[, Rank := seq_len(.N)]
      
      # Calculate effective number of founders before subsetting
      fe2_sum <- sum(founder_dt$Contrib^2, na.rm = TRUE)
      n_founders_total <- nrow(founder_dt)
      
      result$founders <- if (nrow(founder_dt) > top) founder_dt[1:top] else founder_dt
    }
  }
  
  # ---- Ancestor Contributions (Boichard's iterative algorithm) ----
  fa2_sum <- NA_real_
  n_ancestors_total <- 0L
  
  if (mode %in% c("ancestor", "both")) {
    message("Calculating ancestor contributions (Boichard's iterative algorithm)...")
    
    # Work with mutable integer parent arrays
    work_sire <- as.integer(sire_num)
    work_dam <- as.integer(dam_num)
    
    # Find all ancestors of cohort via BFS on integer indices
    is_ancestor <- logical(n)
    visited <- logical(n)
    queue <- integer(n)
    head <- 1L
    tail <- length(cohort_pos)
    if (tail > 0) {
      queue[seq_len(tail)] <- as.integer(cohort_pos)
      visited[queue[seq_len(tail)]] <- TRUE
    }
    while (head <= tail) {
      pos <- queue[head]
      head <- head + 1L
      s <- sire_num[pos]
      d <- dam_num[pos]
      if (s > 0 && !visited[s]) {
        visited[s] <- TRUE
        is_ancestor[s] <- TRUE
        tail <- tail + 1L
        queue[tail] <- s
      }
      if (d > 0 && !visited[d]) {
        visited[d] <- TRUE
        is_ancestor[d] <- TRUE
        tail <- tail + 1L
        queue[tail] <- d
      }
    }
    
    candidate_idx <- which(is_ancestor)
    n_total_anc <- length(candidate_idx)
    
    if (n_total_anc == 0) {
      result$ancestors <- data.table(Ind = character(0), Contrib = numeric(0),
                                      CumContrib = numeric(0), Rank = integer(0))
    } else {
      selected_idx <- integer(n_total_anc)
      marginal_contribs <- numeric(n_total_anc)
      n_selected <- 0L
      active_mask <- rep(TRUE, n_total_anc)
      
      repeat {
        # Backward pass with current (modified) parent links
        contrib <- backward_pass(work_sire, work_dam, cohort_pos, n, n_cohort)
        
        # Find best among remaining candidates
        cand_contribs <- contrib[candidate_idx]
        cand_contribs[!active_mask] <- 0
        
        best_local <- which.max(cand_contribs)
        best_val <- cand_contribs[best_local]
        
        if (best_val <= 1e-10) break
        
        n_selected <- n_selected + 1L
        best_pos <- candidate_idx[best_local]
        selected_idx[n_selected] <- best_pos
        marginal_contribs[n_selected] <- best_val
        
        # "Account for" this ancestor: sever its parental links
        work_sire[best_pos] <- 0L
        work_dam[best_pos] <- 0L
        active_mask[best_local] <- FALSE
      }
      
      if (n_selected == 0) {
        ancestor_dt_full <- data.table(
          Ind = character(0),
          Contrib = numeric(0),
          CumContrib = numeric(0),
          Rank = integer(0)
        )
      } else {
        ancestor_dt_full <- data.table(
          Ind = ind_ids[selected_idx[seq_len(n_selected)]],
          Contrib = marginal_contribs[seq_len(n_selected)]
        )
        ancestor_dt_full[, CumContrib := cumsum(Contrib)]
        ancestor_dt_full[, Rank := seq_len(.N)]
      }
      
      # Calculate effective number of ancestors before subsetting
      fa2_sum <- sum(ancestor_dt_full$Contrib^2, na.rm = TRUE)
      n_ancestors_total <- nrow(ancestor_dt_full)
      
      result$ancestors <- if (nrow(ancestor_dt_full) > top) ancestor_dt_full[1:top] else ancestor_dt_full
    }
  }
  
  # ---- Summary Statistics ----
  summary_list <- list(n_cohort = n_cohort)
  
  if (!is.null(result$founders)) {
    summary_list$n_founders_total <- n_founders_total
    summary_list$n_founders_reported <- nrow(result$founders)
    summary_list$Ne_f <- if (!is.na(fe2_sum) && fe2_sum > 0) 1 / fe2_sum else NA_real_
  }
  
  if (!is.null(result$ancestors)) {
    summary_list$n_ancestors_total <- n_ancestors_total
    summary_list$n_ancestors_reported <- nrow(result$ancestors)
    summary_list$Ne_a <- if (!is.na(fa2_sum) && fa2_sum > 0) 1 / fa2_sum else NA_real_
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
#' @param labelvar Character. The name of the column containing founder labels 
#'   (e.g., "Breed", "Origin").
#' @param labels Character vector. Specific labels to calculate proportions for. 
#'   If NULL, all unique labels in \code{labelvar} for founders are used.
#'
#' @return A \code{data.table} with columns for \code{Ind} and each tracked label.
#'
#' @export
pedancestry <- function(ped, labelvar, labels = NULL) {
  if (!inherits(ped, "tidyped")) stop("ped must be a tidyped object")
  if (!labelvar %in% names(ped)) stop(sprintf("Column '%s' not found.", labelvar))
  
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
  if (is.null(labels)) {
    labels <- unique(ped_work[[labelvar]][founders_idx])
    labels <- labels[!is.na(labels) & labels != ""]
  }
  
  if (length(labels) == 0) {
    stop("No valid labels found in founders. Please check 'labelvar'.")
  }

  n <- nrow(ped_work)
  res_mat <- matrix(0, nrow = n, ncol = length(labels))
  colnames(res_mat) <- labels
  
  # Initialize founders with their respective labels
  for (j in seq_along(labels)) {
    lab <- labels[j]
    idx <- which(ped_work[[labelvar]] == lab & is.na(ped_work$Sire) & is.na(ped_work$Dam))
    res_mat[idx, j] <- 1.0
  }
  
  sire_nums <- ped_work$SireNum
  dam_nums <- ped_work$DamNum
  
  # Map global IndNum to matrix row index in the sorted set
  indnum_to_row <- integer(max(ped_work$IndNum, na.rm = TRUE))
  indnum_to_row[ped_work$IndNum] <- seq_len(n)
  
  # Forward pass via Rcpp
  res_mat <- cpp_calculate_ancestry(
    sire_nums, 
    dam_nums, 
    res_mat, 
    indnum_to_row
  )
  
  ans <- as.data.table(res_mat)
  colnames(ans) <- labels
  ans[, Ind := ped_work$Ind]
  setcolorder(ans, "Ind")
  return(ans)
}

#' Summarize Inbreeding Levels
#'
#' Classifies individuals into inbreeding levels based on their inbreeding 
#' coefficients (F) according to standard academic thresholds.
#'
#' @param ped A \code{tidyped} object.
#'
#' @return A \code{data.table} with counts and percentages for each F-group.
#'
#' @export
pedinbreed_class <- function(ped) {
  if (!inherits(ped, "tidyped")) stop("ped must be a tidyped object")
  
  if (!"f" %in% names(ped)) {
    message("Calculating inbreeding coefficients...")
    ped <- inbreed(ped)
  }
  
  breaks <- c(-Inf, 0.000001, 0.0625, 0.125, 0.25, Inf)
  labels <- c("F = 0", "0 < F <= 0.0625", "0.0625 < F <= 0.125", "0.125 < F <= 0.25", "F > 0.25")
  
  # Classify
  f_classes <- cut(ped$f, breaks = breaks, labels = labels, include.lowest = TRUE)
  
  dt <- data.table(F_Class = factor(f_classes, levels = labels))
  summary_dt <- dt[, .(Count = .N), by = F_Class]
  
  # Ensure all categories appear
  all_cats <- data.table(F_Class = factor(labels, levels = labels))
  summary_dt <- merge(all_cats, summary_dt, by = "F_Class", all.x = TRUE)
  summary_dt[is.na(Count), Count := 0]
  
  summary_dt[, Percentage := (Count / sum(Count)) * 100]
  setorder(summary_dt, F_Class)
  
  return(summary_dt)
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
#' $F_{ij}$. Journal of Heredity, 87(4), 337-339.
#' 
#' Meuwissen, T. H., & Luo, Z. (1992). Computing inbreeding coefficients in 
#' large populations. Genetics Selection Evolution, 24(4), 305-313.
#'
#' @export
pedpartial <- function(ped, ancestors = NULL, top = 20) {
  if (!inherits(ped, "tidyped")) stop("ped must be a tidyped object")
  
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
  
  # Call C++ partial inbreeding engine
  pf_mat <- cpp_calculate_partial_inbreeding(
    ped_work$SireNum, 
    ped_work$DamNum, 
    res_f$dii, 
    as.integer(anc_nums)
  )
  
  result <- as.data.table(pf_mat)
  colnames(result) <- anc_ids
  result[, Ind := ped_work$Ind]
  setcolorder(result, "Ind")
  
  return(result)
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
  cat(sprintf("Reference cohort size: %d\n", s$n_cohort))
  
  if (!is.null(x$founders)) {
    cat(sprintf("\nFounders: %d (reported top %d)\n", s$n_founders_total, s$n_founders_reported))
    cat(sprintf("Effective number of founders (Ne_f): %.2f\n", s$Ne_f))
    cat("\nTop 10 Founder Contributions:\n")
    print(head(x$founders, 10))
  }
  
  if (!is.null(x$ancestors)) {
    cat(sprintf("\nAncestors: %d (reported top %d)\n", s$n_ancestors_total, s$n_ancestors_reported))
    cat(sprintf("Effective number of ancestors (Ne_a): %.2f\n", s$Ne_a))
    cat("\nTop 10 Ancestor Contributions:\n")
    print(head(x$ancestors, 10))
  }
  
  invisible(x)
}
