#' Split Pedigree into Disconnected Groups
#'
#' @description
#' Detects and splits a tidyped object into disconnected groups (connected components).
#' Uses igraph to efficiently find groups of individuals that have no genetic
#' relationships with each other. Isolated individuals (Gen = 0, those with no
#' parents and no offspring) are excluded from group splitting and stored separately.
#'
#' @param ped A tidyped object created by \code{\link{tidyped}}.
#'
#' @return A list of class "splitped" containing:
#'   \item{GP1, GP2, ...}{tidyped objects for each disconnected group (with at least
#'     2 individuals), with renumbered IndNum, SireNum, DamNum}
#'
#'   The returned object has the following attributes:
#'   \item{n_groups}{Number of disconnected groups found (excluding isolated individuals)}
#'   \item{sizes}{Named vector of group sizes}
#'   \item{total}{Total number of individuals in groups (excluding isolated)}
#'   \item{isolated}{Character vector of isolated individual IDs (Gen = 0)}
#'   \item{n_isolated}{Number of isolated individuals}
#'
#' @details
#' This function identifies connected components in the pedigree graph where
#' edges represent parent-offspring relationships. Two individuals are in the
#' same group if they share any ancestry (direct or indirect).
#'
#' Isolated individuals (Gen = 0 in tidyped output) are those who:
#' \itemize{
#'   \item Have no known parents (Sire and Dam are both NA)
#'   \item Are not parents of any other individual in the pedigree
#' }
#' These isolated individuals are excluded from splitting and stored in the
#' \code{isolated} attribute. Each resulting group contains at least 2 individuals
#' (at least one parent-offspring relationship).
#'
#' The function always returns a list, even if there is only one group (i.e.,
#' the pedigree is fully connected). Groups are sorted by size in descending
#' order.
#'
#' Each group in the result is a valid tidyped object with:
#' \itemize{
#'   \item Renumbered IndNum (1 to n for each group)
#'   \item Updated SireNum and DamNum referencing the new IndNum
#'   \item Recalculated Gen (generation) based on the group's structure
#' }
#'
#' @examples
#' # Load example data
#' library(visPedigree)
#' data(small_ped)
#'
#' # First tidy the pedigree
#' tped <- tidyped(small_ped)
#'
#' # Split into groups
#' result <- splitped(tped)
#' print(result)
#'
#' # Access individual groups (each is a tidyped object)
#' result$GP1
#'
#' # Check isolated individuals
#' attr(result, "isolated")
#'
#' @seealso \code{\link{tidyped}} for pedigree tidying
#'
#' @export
splitped <- function(ped) {
  # Input validation: ensure tidyped class (auto-recover if possible)
  ped <- ensure_tidyped(ped)

  # Copy to avoid modifying original
  ped <- data.table::copy(ped)

  # Separate isolated individuals (Gen = 0)
  # These have no parents AND are not parents of anyone
  isolated_idx <- ped$Gen == 0L
  isolated_ids <- ped$Ind[isolated_idx]
  n_isolated <- sum(isolated_idx)

  # Filter to non-isolated individuals only
  ped_valid <- ped[!isolated_idx, ]
  n <- nrow(ped_valid)

  # Handle edge case: no valid individuals (all isolated)
  if (n == 0) {
    result <- list()
    class(result) <- c("splitped", "list")
    attr(result, "n_groups") <- 0L
    attr(result, "sizes") <- integer(0)
    attr(result, "total") <- 0L
    attr(result, "isolated") <- isolated_ids
    attr(result, "n_isolated") <- n_isolated
    return(result)
  }

  # Handle edge case: single non-isolated individual
  # (This shouldn't happen in practice since a valid pedigree needs at least 2)
  if (n == 1) {
    group_ped <- retidy_subset(ped_valid)
    result <- list(GP1 = group_ped)
    class(result) <- c("splitped", "list")
    attr(result, "n_groups") <- 1L
    attr(result, "sizes") <- c(GP1 = 1L)
    attr(result, "total") <- 1L
    attr(result, "isolated") <- isolated_ids
    attr(result, "n_isolated") <- n_isolated
    return(result)
  }

  # Find connected components using numeric indices
  group_ids <- find_components(ped_valid)

  # Build result from group assignments
  build_splitped_result(ped_valid, group_ids, n, isolated_ids, n_isolated)
}

# Find connected components using numeric indices
find_components <- function(ped) {
  n <- nrow(ped)
  ind_nums <- ped$IndNum
  sire_nums <- ped$SireNum
  dam_nums <- ped$DamNum

  # Build edges using numeric indices directly
  # In tidyped, SireNum=0 and DamNum=0 means unknown parent
  # All parents are guaranteed to be in the pedigree

  edges_from <- integer(0)
  edges_to <- integer(0)

  # Add sire-offspring edges
  valid_sire <- sire_nums > 0
  if (any(valid_sire)) {
    edges_from <- c(edges_from, ind_nums[valid_sire])
    edges_to <- c(edges_to, sire_nums[valid_sire])
  }

  # Add dam-offspring edges
  valid_dam <- dam_nums > 0
  if (any(valid_dam)) {
    edges_from <- c(edges_from, ind_nums[valid_dam])
    edges_to <- c(edges_to, dam_nums[valid_dam])
  }

  # Create graph and find components
  if (length(edges_from) == 0) {
    # No edges: each individual is its own group
    return(seq_len(n))
  }

  # Build edge matrix
  edge_matrix <- cbind(edges_from, edges_to)

  # Create graph with numeric vertices
  g <- igraph::graph_from_edgelist(edge_matrix, directed = FALSE)

  # Add isolated vertices if any
  max_vertex <- max(ind_nums)
  if (igraph::vcount(g) < max_vertex) {
    g <- igraph::add_vertices(g, max_vertex - igraph::vcount(g))
  }

  # Find connected components
  components <- igraph::components(g)

  # Map back to original order
  components$membership[ind_nums]
}

# Re-tidy a subset to be a valid standalone tidyped
# Uses tidyped() to properly recalculate Gen, IndNum, SireNum, DamNum
retidy_subset <- function(ped) {
  if (nrow(ped) == 0) return(ped)

  # Extract only the basic pedigree columns (Ind, Sire, Dam, Sex)
  # and call tidyped() to properly re-process the subset
  basic_ped <- data.table::data.table(
    Ind = ped$Ind,
    Sire = ped$Sire,
    Dam = ped$Dam
  )

  # Add Sex column if present
  if ("Sex" %in% names(ped)) {
    basic_ped$Sex <- ped$Sex
  }

  # Check if all parents are missing (founders only group)
  # tidyped() does not accept pedigrees where all parents are NA
  all_founders <- all(is.na(basic_ped$Sire) & is.na(basic_ped$Dam))

  if (all_founders) {
    # Build a minimal tidyped manually for founders-only groups
    n <- nrow(basic_ped)
    result <- data.table::data.table(
      Ind = basic_ped$Ind,
      Sire = rep(NA_character_, n),
      Dam = rep(NA_character_, n),
      Gen = rep(1L, n),
      Sex = if ("Sex" %in% names(basic_ped)) basic_ped$Sex else rep(NA_character_, n),
      IndNum = seq_len(n),
      SireNum = rep(0L, n),
      DamNum = rep(0L, n)
    )
    result <- new_tidyped(result)
    meta <- attr(ped, "ped_meta")
    if (!is.null(meta)) data.table::setattr(result, "ped_meta", meta)
    return(result)
  }

  # Propagate selfing attribute if present in the original pedigree
  meta <- attr(ped, "ped_meta")
  selfing_val <- if (!is.null(meta)) isTRUE(meta$selfing) else FALSE
  
  # Call tidyped to properly recalculate Gen, IndNum, SireNum, DamNum
  tidyped(basic_ped, selfing = selfing_val)
}

# Build result from group membership vector
build_splitped_result <- function(ped, group_ids, n, isolated_ids, n_isolated) {
  # Get unique groups and their sizes
  group_table <- table(group_ids)
  group_order <- order(group_table, decreasing = TRUE)
  unique_groups <- as.integer(names(group_table)[group_order])

  n_groups <- length(unique_groups)

  # Create result list
  result <- vector("list", n_groups)
  sizes <- integer(n_groups)

  for (i in seq_len(n_groups)) {
    grp <- unique_groups[i]
    idx <- which(group_ids == grp)
    group_ped <- ped[idx, ]

    # Re-tidy the subset using tidyped() to properly recalculate Gen, IndNum, etc.
    group_ped <- retidy_subset(group_ped)

    result[[i]] <- group_ped
    sizes[i] <- length(idx)
  }

  # Name the groups
  names(result) <- paste0("GP", seq_len(n_groups))
  names(sizes) <- names(result)

  # Set class and attributes
  class(result) <- c("splitped", "list")
  attr(result, "n_groups") <- n_groups
  attr(result, "sizes") <- sizes
  attr(result, "total") <- n
  attr(result, "isolated") <- isolated_ids
  attr(result, "n_isolated") <- n_isolated

  result
}

#' @export
print.splitped <- function(x, ...) {
  n_groups <- attr(x, "n_groups")
  sizes <- attr(x, "sizes")
  total <- attr(x, "total")
  n_isolated <- attr(x, "n_isolated")

  cat("Pedigree Split Result\n")
  cat("======================\n")
  cat("Total individuals in groups:", total, "\n")
  cat("Isolated individuals (Gen=0):", n_isolated, "\n")
  cat("Number of groups: ", n_groups, "\n")

  if (n_groups > 0) {
    cat("\nGroup sizes:\n")
    for (i in seq_len(min(n_groups, 10))) {
      cat(sprintf("  %s: %d individuals (%.1f%%)\n",
                  names(sizes)[i], sizes[i], 100 * sizes[i] / total))
    }
    if (n_groups > 10) {
      cat(sprintf("  ... and %d more groups\n", n_groups - 10))
    }
  }

  invisible(x)
}

#' @export
summary.splitped <- function(object, ...) {
  n_groups <- attr(object, "n_groups")
  sizes <- attr(object, "sizes")
  total <- attr(object, "total")
  n_isolated <- attr(object, "n_isolated")

  cat("Summary of Pedigree Split\n")
  cat("=========================\n")
  cat("Total individuals in groups:", total, "\n")
  cat("Isolated individuals (Gen=0):", n_isolated, "\n")
  cat("Grand total:", total + n_isolated, "\n")
  cat("Number of groups: ", n_groups, "\n")

  if (n_groups > 0) {
    cat("\nSize statistics:\n")
    cat("  Min:    ", min(sizes), "\n")
    cat("  Max:    ", max(sizes), "\n")
    cat("  Mean:   ", round(mean(sizes), 1), "\n")
    cat("  Median: ", median(sizes), "\n")

    if (n_groups > 1) {
      cat("\nConnectivity: Pedigree contains disconnected groups\n")
    } else {
      cat("\nConnectivity: Pedigree is fully connected\n")
    }
  }

  invisible(object)
}
