#' Test if an object is a tidyped
#'
#' @param x An object to test.
#' @return Logical scalar.
#' @export
is_tidyped <- function(x) {
  inherits(x, "tidyped")
}

#' Build a unified metadata list for a tidyped object
#'
#' @param selfing Logical: whether selfing/monoecious mode was used.
#' @param bisexual_parents Character vector of IDs that appear as both sire and dam.
#' @param genmethod Character: generation assignment method ("top" or "bottom").
#' @return A named list of pedigree metadata.
#' @keywords internal
build_ped_meta <- function(selfing = FALSE,
                           bisexual_parents = character(0),
                           genmethod = "top") {
  list(
    selfing           = selfing,
    bisexual_parents  = bisexual_parents,
    genmethod         = genmethod
  )
}

#' Access pedigree metadata from a tidyped object
#'
#' @param x A tidyped object.
#' @return The \code{ped_meta} list, or \code{NULL} if not set.
#' @export
pedmeta <- function(x) {
  if (!is_tidyped(x)) stop("x is not a tidyped object.", call. = FALSE)
  attr(x, "ped_meta")
}

#' Check whether a tidyped object contains inbreeding coefficients
#'
#' @param x A tidyped object.
#' @return Logical scalar.
#' @export
has_inbreeding <- function(x) {
  is_tidyped(x) && "f" %in% names(x)
}

#' Check whether a tidyped object contains candidate flags
#'
#' @param x A tidyped object.
#' @return Logical scalar.
#' @export
has_candidates <- function(x) {
  is_tidyped(x) && "Cand" %in% names(x)
}

#' Internal constructor for tidyped class
#' @param x A data.table object
#' @return A tidyped object
#' @keywords internal
new_tidyped <- function(x) {
  stopifnot(inherits(x, "data.table"))
  if (!is_tidyped(x)) {
    data.table::setattr(x, "class", c("tidyped", class(x)))
  }
  # x[] clears data.table's invisible flag set by := and set*() operations.
  # [.tidyped uses setattr (in-place), so no copy is created.
  x[]
}

#' Check whether all referenced parents are present
#'
#' @param x A pedigree-like object.
#' @return Logical scalar.
#' @keywords internal
is_complete_pedigree <- function(x) {
  if (!is.data.frame(x)) return(FALSE)

  core_cols <- c("Ind", "Sire", "Dam")
  if (!all(core_cols %in% names(x))) return(FALSE)

  sires_ok <- is.na(x$Sire) | x$Sire %in% x$Ind
  dams_ok  <- is.na(x$Dam)  | x$Dam %in% x$Ind

  all(sires_ok) && all(dams_ok)
}

#' Error on row-truncated pedigree subsets
#'
#' @param ped A pedigree-like object.
#' @param fun Character scalar. Calling function name for the error message.
#' @return The input object, invisibly.
#' @keywords internal
require_complete_pedigree <- function(ped, fun) {
  if (is_complete_pedigree(ped)) return(invisible(ped))

  stop(
    fun, " requires a structurally complete pedigree. ",
    "This input appears to be a row-truncated subset with missing parent records.\n",
    "Compute on the full pedigree first, or extract a valid sub-pedigree with ",
    "`tidyped(tp, cand = ids, trace = \"up\")`.",
    call. = FALSE
  )
}

#' Internal helper to ensure ped is a tidyped object
#'
#' If the object has lost its \code{tidyped} class (e.g., after \code{merge()},
#' \code{rbind()}, or dplyr operations) but still contains the required columns,
#' the class is automatically restored with an informational message. Otherwise,
#' an error is raised guiding the user to call \code{tidyped()} or
#' \code{as_tidyped()}.
#'
#' @param ped An object expected to be a tidyped.
#' @return A valid tidyped object.
#' @keywords internal
ensure_tidyped <- function(ped) {
  if (is_tidyped(ped)) return(ped)

  if (!is.data.frame(ped)) {
    stop("'ped' must be a tidyped object. Run `tp <- tidyped(ped)` first.",
         call. = FALSE)
  }

  core <- c("Ind", "Sire", "Dam", "Sex", "Gen", "IndNum", "SireNum", "DamNum")
  missing <- setdiff(core, names(ped))

  if (length(missing) > 0) {
    stop("'ped' must be a tidyped object. Run `tp <- tidyped(ped)` first.",
         call. = FALSE)
  }

  if (!inherits(ped, "data.table")) {
    ped <- data.table::as.data.table(ped)
  }

  # Validate that IndNum indices are still consistent with Ind/Sire/Dam.
  # After rbind() or row subsetting, IndNum may be duplicated or no longer
  # match the row position, which would cause silent errors in Rcpp routines.
  if (anyDuplicated(ped$IndNum) > 0L ||
      !identical(ped$IndNum, seq_len(nrow(ped)))) {
    # Rebuild numeric indices from current row order
    ped[, IndNum := .I]
    ped[, SireNum := match(Sire, Ind, nomatch = 0L)]
    ped[, DamNum  := match(Dam, Ind, nomatch = 0L)]
  }

  # Structure is intact, only the class label was dropped
  message(
    "Note: 'ped' lost its tidyped class ",
    "(common after merge/rbind/dplyr). Restoring automatically."
  )
  new_tidyped(ped)
}

#' Internal helper to ensure ped is a complete tidyped object
#'
#' Like \\code{ensure_tidyped()}, but also rejects row-truncated pedigree
#' subsets whose referenced parents are no longer present.
#'
#' @param ped An object expected to be a complete tidyped pedigree.
#' @param fun Character scalar. Calling function name for the error message.
#' @return A valid, structurally complete tidyped object.
#' @keywords internal
ensure_complete_tidyped <- function(ped, fun) {
  if (is_tidyped(ped)) {
    require_complete_pedigree(ped, fun)
    return(ped)
  }

  if (!is.data.frame(ped)) {
    stop("'ped' must be a tidyped object. Run `tp <- tidyped(ped)` first.",
         call. = FALSE)
  }

  core <- c("Ind", "Sire", "Dam", "Sex", "Gen", "IndNum", "SireNum", "DamNum")
  missing <- setdiff(core, names(ped))

  if (length(missing) > 0) {
    stop("'ped' must be a tidyped object. Run `tp <- tidyped(ped)` first.",
         call. = FALSE)
  }

  if (!inherits(ped, "data.table")) {
    ped <- data.table::as.data.table(ped)
  }

  require_complete_pedigree(ped, fun)

  if (anyDuplicated(ped$IndNum) > 0L ||
      !identical(ped$IndNum, seq_len(nrow(ped)))) {
    ped[, IndNum := .I]
    ped[, SireNum := match(Sire, Ind, nomatch = 0L)]
    ped[, DamNum  := match(Dam, Ind, nomatch = 0L)]
  }

  message(
    "Note: 'ped' lost its tidyped class ",
    "(common after merge/rbind/dplyr). Restoring automatically."
  )
  new_tidyped(ped)
}

#' Restore the tidyped class to a manipulated pedigree
#'
#' Rapidly restores the \code{tidyped} class to a \code{data.table} or
#' \code{data.frame} that was previously processed by \code{\link{tidyped}()}
#' but lost its class attributes due to data manipulation.
#'
#' This is a lightweight operation that only checks for the required columns
#' and re-attaches the class---it does \strong{not} re-run the full pedigree
#' sorting, generation inference, or loop detection.
#'
#' @param x A \code{data.table} or \code{data.frame} that was previously a
#'   tidyped object. It must still contain the core columns: \code{Ind},
#'   \code{Sire}, \code{Dam}, \code{Sex}, \code{Gen}, \code{IndNum},
#'   \code{SireNum}, \code{DamNum}.
#'
#' @return A \code{tidyped} object.
#'
#' @details
#' This helper is intended for objects that still contain the core pedigree
#' columns and numeric indices, but no longer inherit from \code{tidyped}.
#' A common reproducible case is \code{rbind()} on two \code{tidyped}
#' fragments, which typically returns a plain \code{data.table}. Converting a
#' \code{tidyped} object to a plain \code{data.frame} and then subsetting it
#' also drops the class.
#'
#' Some operations, such as \code{merge()} or certain dplyr workflows, may or
#' may not preserve the \code{tidyped} class depending on the versions of
#' \pkg{data.table}, \pkg{dplyr}, and the exact method dispatch path used in
#' the current R session. Therefore, \code{as_tidyped()} should be viewed as a
#' safe recovery helper rather than something only needed after one specific
#' verb.
#'
#' Typical class-loss scenarios include:
#' \itemize{
#'   \item \code{rbind(tped1, tped2)} --- often returns plain \code{data.table}
#'   \item \code{as.data.frame(tped)[rows, ]} --- returns plain \code{data.frame}
#'   \item manual class removal or serialization / import workflows
#' }
#' After such operations, downstream analysis functions (e.g.,
#' \code{\link{pedstats}}, \code{\link{pedne}}) will either error or
#' automatically restore the class. You can also call \code{as_tidyped()}
#' explicitly to restore the class yourself.
#'
#' @examples
#' library(visPedigree)
#' tp <- tidyped(simple_ped)
#' class(tp)
#' # [1] "tidyped"    "data.table" "data.frame"
#'
#' # Simulate class loss via rbind()
#' tp2 <- rbind(tp[1:5], tp[6:10])
#' class(tp2)
#' # [1] "data.table" "data.frame"
#'
#' # Restore the class
#' tp3 <- as_tidyped(tp2)
#' class(tp3)
#' # [1] "tidyped"    "data.table" "data.frame"
#'
#' # It can also restore from a plain data.frame if core columns are intact
#' tp_df <- as.data.frame(tp)
#' tp4 <- tp_df[tp_df$Gen > 1, ]
#' class(tp4)
#' # [1] "data.frame"
#'
#' tp5 <- as_tidyped(tp4)
#' class(tp5)
#' # [1] "tidyped"    "data.table" "data.frame"
#'
#' @seealso \code{\link{tidyped}}, \code{\link{new_tidyped}}
#' @export
as_tidyped <- function(x) {
  if (is_tidyped(x)) return(x)

  if (!is.data.frame(x)) {
    stop("Cannot coerce to tidyped: input is not a data.frame or data.table.",
         call. = FALSE)
  }

  if (!inherits(x, "data.table")) {
    x <- data.table::as.data.table(x)
  }

  core <- c("Ind", "Sire", "Dam", "Sex", "Gen", "IndNum", "SireNum", "DamNum")
  missing <- setdiff(core, names(x))

  if (length(missing) > 0) {
    stop(
      "Cannot restore tidyped class. Missing columns: ",
      paste(missing, collapse = ", "),
      ". Run `tidyped()` on raw data instead.",
      call. = FALSE
    )
  }

  # Rebuild numeric indices if they are inconsistent (e.g. after rbind/subset)
  if (anyDuplicated(x$IndNum) > 0L ||
      !identical(x$IndNum, seq_len(nrow(x)))) {
    x[, IndNum := .I]
    x[, SireNum := match(Sire, Ind, nomatch = 0L)]
    x[, DamNum  := match(Dam, Ind, nomatch = 0L)]
  }

  new_tidyped(x)
}

#' Internal validator for tidyped class
#'
#' Validates a tidyped object. If the object has lost its class but retains the
#' required columns, it is automatically restored via \code{ensure_tidyped()}.
#' Fatal structural problems (e.g., missing core columns) raise an error.
#'
#' @param x A tidyped object
#' @return The (possibly restored) object if valid, otherwise an error
#' @keywords internal
validate_tidyped <- function(x) {
  # Attempt auto-recovery if class was lost

  if (!is_tidyped(x)) {
    x <- tryCatch(
      suppressMessages(ensure_tidyped(x)),
      error = function(e) {
        stop("Object must be of class 'tidyped'. ", conditionMessage(e),
             call. = FALSE)
      }
    )
  }

  needed <- c("Ind", "Sire", "Dam", "Sex")
  missing_cols <- setdiff(needed, names(x))
  
  if (length(missing_cols) > 0) {
    stop(
      "The pedigree object is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      ". Please ensure it has been processed by tidyped().",
      call. = FALSE
    )
  }
  
  x
}

#' Subset a tidyped object
#'
#' Intercepts \code{data.table}'s \code{[} method for \code{tidyped} objects.
#' After subsetting, the method checks whether the result is still a valid
#' pedigree (all referenced parents still present). If so, \code{IndNum},
#' \code{SireNum}, and \code{DamNum} are rebuilt and the \code{tidyped} class
#' is preserved. If the pedigree becomes structurally incomplete (missing parent
#' records), the result is degraded to a plain \code{data.table} with a warning.
#' Column-only selections (missing core columns) also return a plain
#' \code{data.table}.
#'
#' @param x A \code{tidyped} object.
#' @param i,j,... Arguments passed to the \code{data.table} \code{[} method.
#' @return A \code{tidyped} object if the result is still a complete pedigree,
#'   otherwise a plain \code{data.table}.
#' @export
`[.tidyped` <- function(x, ...) {
  # Save metadata before dispatching to data.table
  meta <- attr(x, "ped_meta")
  cl <- class(x)

  sc <- sys.call()

  # Detect := (modify-in-place) operations.
  # := may appear in any argument position, e.g.:
  #   x[, col := val]              -> sc[[3]] is `:=`(col, val)
  #   x[cond, col := val]          -> sc[[4]] is `:=`(col, val)
  #   x[, `:=`(col1=v1, col2=v2)] -> sc[[3]] is `:=`(...)
  has_assign <- FALSE
  for (k in seq_along(sc)[-c(1L, 2L)]) {
    is_assign_k <- tryCatch({
      arg <- sc[[k]]
      is.call(arg) && identical(arg[[1L]], quote(`:=`))
    }, error = function(e) FALSE)
    if (isTRUE(is_assign_k)) {
      has_assign <- TRUE
      break
    }
  }

  # Use data.table::setattr() to strip the class IN PLACE (no copy!).
  # This is critical: class(x) <- ... triggers copy-on-modify, which breaks

  # := operations because the caller's variable still points to the original.
  data.table::setattr(x, "class", setdiff(cl, "tidyped"))
  # Ensure class is always restored, even if an error occurs
  on.exit(data.table::setattr(x, "class", cl), add = TRUE)

  # Re-dispatch the original call with proper NSE support.
  # NextMethod() breaks data.table's substitute() mechanism, so we
  # reconstruct the call and evaluate in a helper environment that
  # chains to the caller's frame (for user-defined variables in i/j).
  sc[[1L]] <- quote(`[`)
  sc[[2L]] <- quote(x)
  env <- list2env(list(x = x), parent = parent.frame())
  result <- eval(sc, envir = env)

  # Restore class on the original object (in-place)
  data.table::setattr(x, "class", cl)

  if (has_assign) {
    # := modifies x by reference. x IS the caller's object (no copy).
    # Class already restored above. Return invisibly.
    return(invisible(x))
  }

  # If result is not a data.table (e.g., single column vector), return as-is
  if (!is.data.table(result)) return(result)

  # Check if core columns are still present (column subset may have removed them)
  core_cols <- c("Ind", "Sire", "Dam")
  if (!all(core_cols %in% names(result))) return(result)

  # Check pedigree completeness: all non-NA parents must still be in Ind
  if (is_complete_pedigree(result)) {
    # Pedigree is still complete -> rebuild IndNum and restore class.
    # Use data.table::set() to avoid recursive [.tidyped dispatch.
    if ("IndNum" %in% names(result)) {
      data.table::set(result, j = "IndNum",  value = seq_len(nrow(result)))
      data.table::set(result, j = "SireNum",
                      value = match(result$Sire, result$Ind, nomatch = 0L))
      data.table::set(result, j = "DamNum",
                      value = match(result$Dam, result$Ind, nomatch = 0L))
    }
    data.table::setattr(result, "ped_meta", meta)
    data.table::setattr(result, "class", cl)
  } else {
    # Pedigree is incomplete -> degrade to plain data.table
    warning(
      "Subsetting removed parent records. ",
      "Result is a plain data.table, not a tidyped.\n",
      "Use tidyped(tp, cand = ids, trace = \"up\") ",
      "to extract a valid sub-pedigree.",
      call. = FALSE
    )
  }

  result
}

#' Print method for tidyped pedigree
#' @param x A tidyped object
#' @param ... Additional arguments passed to the data.table print method
#' @return The input object, invisibly.
#' @export
print.tidyped <- function(x, ...) {
  cat("Tidy Pedigree Object\n")
  # Temporarily remove tidyped class to use data.table's print method
  cl <- class(x)
  class(x) <- setdiff(cl, "tidyped")
  print(x, ...)
  # Restore the class
  class(x) <- cl
  invisible(x)
}

#' Plot a tidy pedigree
#'
#' @param x A \code{tidyped} object.
#' @param ... Additional arguments passed to \code{\link{visped}}.
#' @return Invisibly returns a list of graph data from \code{\link{visped}} (node/edge
#'   data and layout components) used to render the pedigree; the primary result is
#'   the plot drawn on the current device.
#' @export
plot.tidyped <- function(x, ...) {
  validate_tidyped(x)
  visped(x, ...)
}

#' Summary method for tidyped objects
#'
#' @param object A tidyped object.
#' @param ... Additional arguments (ignored).
#' @return A summary.tidyped object (list) containing:
#' \itemize{
#'   \item \code{n_ind}: Total number of individuals.
#'   \item \code{n_male}, \code{n_female}, \code{n_unknown_sex}: Sex composition counts.
#'   \item \code{n_founders}: Number of individuals with no known parents.
#'   \item \code{n_both_parents}: Count of individuals with complete parentage.
#'   \item \code{max_gen}, \code{gen_dist}: (Optional) Maximum generation and its distribution.
#'   \item \code{n_families}, \code{family_sizes}, \code{top_families}: (Optional) Family statistics.
#'   \item \code{f_stats}, \code{n_inbred}: (Optional) Inbreeding coefficient statistics.
#'   \item \code{n_cand}, \code{cand_f_stats}: (Optional) Candidate-specific statistics.
#' }
#' @export
summary.tidyped <- function(object, ...) {
  x <- validate_tidyped(object)
  res <- list()
  
  # Basic counts
  res$n_ind <- nrow(x)
  res$n_male <- sum(x$Sex == "male", na.rm = TRUE)
  res$n_female <- sum(x$Sex == "female", na.rm = TRUE)
  res$n_monoecious <- sum(x$Sex == "monoecious", na.rm = TRUE)
  res$n_unknown_sex <- sum(is.na(x$Sex))
  
  # Founders: Individuals with both parents unknown (NA)
  res$n_founders <- sum(is.na(x$Sire) & is.na(x$Dam))
  
  # Pedigree completeness
  res$n_both_parents <- sum(!is.na(x$Sire) & !is.na(x$Dam))
  res$n_sire_only <- sum(!is.na(x$Sire) & is.na(x$Dam))
  res$n_dam_only <- sum(is.na(x$Sire) & !is.na(x$Dam))
  
  # Isolated: Gen == 0
  if ("Gen" %in% names(x)) {
    res$n_isolated <- sum(x$Gen == 0)
    res$max_gen <- max(x$Gen, na.rm = TRUE)
    
    # Generation distribution
    gen_table <- table(x$Gen)
    res$gen_dist <- as.data.frame(gen_table)
    names(res$gen_dist) <- c("Generation", "Count")
  } else {
    res$n_isolated <- NA
    res$max_gen <- NA
    res$gen_dist <- NULL
  }
  
  # Candidates
  if ("Cand" %in% names(x)) {
    res$n_cand <- sum(x$Cand, na.rm = TRUE)
  } else {
    res$n_cand <- NA
  }

  # Offspring statistics
  # Count offspring for each individual
  if (!is.null(x$Sire) && !is.null(x$Dam)) {
    sire_counts <- table(x$Sire[!is.na(x$Sire)])
    dam_counts <- table(x$Dam[!is.na(x$Dam)])
    
    # Individuals with offspring (as sire)
    res$n_sires <- length(sire_counts)
    res$max_offspring_sire <- if(length(sire_counts) > 0) max(sire_counts) else 0
    res$mean_offspring_sire <- if(length(sire_counts) > 0) mean(sire_counts) else 0
    
    # Individuals with offspring (as dam)
    res$n_dams <- length(dam_counts)
    res$max_offspring_dam <- if(length(dam_counts) > 0) max(dam_counts) else 0
    res$mean_offspring_dam <- if(length(dam_counts) > 0) mean(dam_counts) else 0
    
    # Combined: any individual with offspring
    parents <- unique(c(names(sire_counts), names(dam_counts)))
    res$n_parents <- length(parents)
  }

  # Family information
  if ("Family" %in% names(x)) {
    # Count non-NA families
    families <- x$Family[!is.na(x$Family)]
    if (length(families) > 0) {
      family_table <- table(families)
      res$n_families <- length(family_table)
      res$family_sizes <- as.numeric(family_table)
      res$max_family_size <- max(res$family_sizes)
      res$mean_family_size <- mean(res$family_sizes)
      
      # Top 5 largest families
      top5 <- head(sort(family_table, decreasing = TRUE), 5)
      res$top_families <- data.frame(
        Family = names(top5),
        Size = as.numeric(top5),
        stringsAsFactors = FALSE
      )
    } else {
      res$n_families <- 0
      res$top_families <- NULL
    }
  }

  # Inbreeding coefficients statistics
  if ("f" %in% names(x)) {
    has_f <- any(!is.na(x$f))
    if (has_f) {
      res$f_stats <- list(
        min = min(x$f, na.rm = TRUE),
        max = max(x$f, na.rm = TRUE),
        mean = mean(x$f, na.rm = TRUE)
      )
    } else {
      res$f_stats <- list(
        min = NA_real_,
        max = NA_real_,
        mean = NA_real_
      )
    }
    
    # Count inbred individuals (f > 0)
    res$n_inbred <- sum(x$f > 0, na.rm = TRUE)

    # Stats for candidates if present
    if (!is.na(res$n_cand) && res$n_cand > 0) {
      cand_f <- x[x$Cand == TRUE, f]
      if (any(!is.na(cand_f))) {
        res$cand_f_stats <- list(
          min = min(cand_f, na.rm = TRUE),
          max = max(cand_f, na.rm = TRUE),
          mean = mean(cand_f, na.rm = TRUE)
        )
        res$n_cand_inbred <- sum(cand_f > 0, na.rm = TRUE)
      }
    }
  }

  class(res) <- "summary.tidyped"
  res
}

#' Print method for summary.tidyped
#'
#' @param x A summary.tidyped object.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @export
print.summary.tidyped <- function(x, ...) {
  cat("Pedigree Summary\n")
  cat("================\n\n")
  
  # Basic information
  cat("Total Individuals: ", x$n_ind, "\n")
  pct_male <- if(x$n_ind > 0) 100 * x$n_male / x$n_ind else 0
  pct_female <- if(x$n_ind > 0) 100 * x$n_female / x$n_ind else 0
  
  cat("  - Males:   ", sprintf("%d (%.1f%%)", x$n_male, pct_male), "\n")
  cat("  - Females: ", sprintf("%d (%.1f%%)", x$n_female, pct_female), "\n")
  
  if (!is.null(x$n_monoecious) && x$n_monoecious > 0) {
    pct_monoecious <- if(x$n_ind > 0) 100 * x$n_monoecious / x$n_ind else 0
    cat("  - Monoecious: ", sprintf("%d (%.1f%%)", x$n_monoecious, pct_monoecious), "\n")
  }
  
  if (x$n_unknown_sex > 0) {
    pct_unknown <- if(x$n_ind > 0) 100 * x$n_unknown_sex / x$n_ind else 0
    cat("  - Unknown: ", sprintf("%d (%.1f%%)", x$n_unknown_sex, pct_unknown), "\n")
  }
  cat("\n")
  
  # Pedigree structure
  cat("Pedigree Structure:\n")
  cat("  - Founders (no parents):  ", x$n_founders, "\n")
  cat("  - Both parents known:     ", x$n_both_parents, "\n")
  if (x$n_sire_only > 0) {
    cat("  - Sire only known:        ", x$n_sire_only, "\n")
  }
  if (x$n_dam_only > 0) {
    cat("  - Dam only known:         ", x$n_dam_only, "\n")
  }
  
  if (!is.na(x$n_isolated) && x$n_isolated > 0) {
    cat("  - Isolated (Gen 0):       ", x$n_isolated, "\n")
  }
  cat("\n")
  
  # Generation information
  if (!is.na(x$max_gen)) {
    cat("Generation:\n")
    cat("  - Maximum: ", x$max_gen, "\n")
    if (!is.null(x$gen_dist) && nrow(x$gen_dist) <= 10) {
      cat("  - Distribution:\n")
      for (i in seq_len(nrow(x$gen_dist))) {
        cat(sprintf("      Gen %s: %d individuals\n", 
                    x$gen_dist[i, 1], x$gen_dist[i, 2]))
      }
    } else if (!is.null(x$gen_dist)) {
      cat(sprintf("  - Distribution: %d generations (use str() for details)\n", 
                  nrow(x$gen_dist)))
    }
    cat("\n")
  }
  
  # Reproduction information
  if (!is.null(x$n_parents)) {
    cat("Reproduction:\n")
    cat("  - Individuals with offspring: ", x$n_parents, "\n")
    cat("  - Sires: ", x$n_sires, 
        sprintf(" (Mean=%.1f, Max=%d offspring)\n", 
                x$mean_offspring_sire, x$max_offspring_sire))
    cat("  - Dams:  ", x$n_dams, 
        sprintf(" (Mean=%.1f, Max=%d offspring)\n", 
                x$mean_offspring_dam, x$max_offspring_dam))
    cat("\n")
  }
  
  # Family information
  if (!is.null(x$n_families) && x$n_families > 0) {
    cat("Full-sibling Families:\n")
    cat("  - Number of families:     ", x$n_families, "\n")
    cat("  - Mean family size:       ", sprintf("%.2f\n", x$mean_family_size))
    cat("  - Maximum family size:    ", x$max_family_size, "\n")
    
    if (!is.null(x$top_families) && nrow(x$top_families) > 0) {
      cat("  - Top families by size:\n")
      for (i in seq_len(nrow(x$top_families))) {
        cat(sprintf("      %s: %d\n", 
                    x$top_families$Family[i], 
                    x$top_families$Size[i]))
      }
    }
    cat("\n")
  }
  
  # Candidates
  if (!is.na(x$n_cand)) {
    cat("Candidates Traced: ", x$n_cand, "\n\n")
  }

  # Inbreeding information
  if (!is.null(x$f_stats)) {
    cat("Inbreeding Coefficients:\n")
    if (all(is.na(unlist(x$f_stats)))) {
      cat("  - (Calculated column present, but all values are NA)\n")
    } else {
      cat("  - All individuals:\n")
      cat(sprintf("      Mean = %.4f, Min = %.4f, Max = %.4f\n",
                  x$f_stats$mean, x$f_stats$min, x$f_stats$max))
      if (!is.null(x$n_inbred)) {
        pct_inbred <- if(x$n_ind > 0) 100 * x$n_inbred / x$n_ind else 0
        cat(sprintf("      Inbred (f > 0): %d (%.1f%%)\n",
                    x$n_inbred, pct_inbred))
      }
      
      if (!is.null(x$cand_f_stats)) {
        cat("  - Candidates:\n")
        cat(sprintf("      Mean = %.4f, Min = %.4f, Max = %.4f\n",
                    x$cand_f_stats$mean, x$cand_f_stats$min, x$cand_f_stats$max))
        if (!is.null(x$n_cand_inbred)) {
          pct_cand_inbred <- if(x$n_cand > 0) 100 * x$n_cand_inbred / x$n_cand else 0
          cat(sprintf("      Inbred (f > 0): %d (%.1f%%)\n",
                      x$n_cand_inbred, pct_cand_inbred))
        }
      }
    }
    cat("\n")
  }
  
  cat("================\n")
  invisible(x)
}
