#' Internal constructor for tidyped class
#' @param x A data.table object
#' @return A tidyped object
#' @keywords internal
new_tidyped <- function(x) {
  stopifnot(inherits(x, "data.table"))
  attr(x, "tidyped") <- TRUE
  if (!inherits(x, "tidyped")) {
    class(x) <- c("tidyped", class(x))
  }
  x[]
}

#' Internal validator for tidyped class
#' @param x A tidyped object
#' @return The object if valid, otherwise an error
#' @keywords internal
validate_tidyped <- function(x) {
  if (!inherits(x, "tidyped")) {
    stop("Object must be of class 'tidyped'", call. = FALSE)
  }
  if (!inherits(x, "data.table")) {
    stop("Object must inherit from 'data.table'", call. = FALSE)
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
  
  if (!isTRUE(attr(x, "tidyped"))) {
    stop(
      "The object must have the 'tidyped' attribute set to TRUE.",
      call. = FALSE
    )
  }
  
  x
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
