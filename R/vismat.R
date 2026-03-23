# ---------------------------------------------------------------------------
# vismat internal thresholds
# ---------------------------------------------------------------------------
# Centralised here so that documentation, code, and tests stay consistent.
# All values are the *inclusive upper bound* for the "lighter" behaviour,
# i.e. the condition that triggers the heavier/simpler fallback is "> MAX".

VISMAT_REORDER_MAX      <- 2000L   # N above this: skip hierarchical clustering
VISMAT_LABEL_MAX        <- 500L    # N above this: hide individual labels
VISMAT_GRID_MAX         <- 100L    # N above this: hide cell grid lines
VISMAT_EXPAND_MAX       <- 5000L   # original N above this: use compact representative view
VISMAT_BLOCK_AGG_MIN_N  <- 10000L  # aggregation switches to block algorithm when
VISMAT_BLOCK_AGG_MIN_G  <- 100L    #   n > _MIN_N  AND  n_grp > _MIN_G

#' Visualize Relationship Matrices
#'
#' @description
#' \code{vismat} provides visualization tools for relationship matrices (A, D, AA),
#' supporting individual-level heatmaps and relationship coefficient histograms.
#' This function is useful for exploring population genetic structure, identifying
#' inbred individuals, and analyzing kinship between families.
#'
#' @param mat A relationship matrix. Can be one of the following types:
#' \itemize{
#'   \item A \code{pedmat} object returned by \code{\link{pedmat}} — including
#'         compact matrices. When \code{by} is specified, group-level means are
#'         computed directly from the compact matrix (no full expansion needed).
#'         Without \code{by}, compact matrices are automatically expanded to
#'         full dimensions before plotting (see Details).
#'   \item A \code{\link{tidyped}} object (automatically calculates additive
#'         relationship matrix A).
#'   \item A named list containing matrices (preferring A, D, AA).
#'   \item A standard \code{matrix} or \code{Matrix} object.
#' }
#' \strong{Note}: Inverse matrices (Ainv, Dinv, AAinv) are not supported for 
#' visualization because their elements do not represent meaningful relationship
#' coefficients.
#'
#' @param ped Optional. A tidied pedigree object (\code{tidyped}), used for 
#' extracting labels or grouping information. Required when using the 
#' \code{by} parameter with a plain matrix input. If \code{mat} is a 
#' \code{pedmat} object, the pedigree is extracted automatically.
#'
#' @param type Character, type of visualization. Supported options:
#' \itemize{
#'   \item \code{"heatmap"}: Relationship matrix heatmap (default). Uses a 
#'         Nature Genetics style color palette (white-orange-red-dark red),
#'         with optional hierarchical clustering and group aggregation.
#'   \item \code{"histogram"}: Distribution histogram of relationship coefficients.
#'         Shows the frequency distribution of lower triangular elements (pairwise kinship).
#' }
#'
#' @param ids Character vector specifying individual IDs to display. Used to 
#' filter and display a submatrix of specific individuals. If \code{NULL} 
#' (default), all individuals are shown.
#'
#' @param reorder Logical. If \code{TRUE} (default), rows and columns are 
#' reordered using hierarchical clustering (Ward.D2 method) to bring closely 
#' related individuals together. Only affects heatmap visualization. 
#' Automatically skipped for large matrices (N > VISMAT_REORDER_MAX, default 2 000)
#' to improve performance.
#' 
#' \strong{Clustering principle}: Based on relationship profile distance (Euclidean
#' distance between rows). Full-sibs have nearly identical relationship profiles
#' with the whole population, so they cluster tightly together and appear as
#' contiguous blocks in the heatmap.
#'
#' @param by Optional. Column name in \code{ped} to group by (e.g., 
#' \code{"Family"}, \code{"Gen"}, \code{"Year"}). When grouping is enabled:
#' \itemize{
#'   \item Individual-level matrix is aggregated to a group-level matrix 
#'         (computing mean relationship coefficients between groups).
#'   \item For \code{"Family"} grouping, founders without family assignment are excluded.
#'   \item For other grouping columns, NA values are assigned to an \code{"Unknown"} group.
#' }
#' Useful for visualizing population structure in large pedigrees.
#'
#' @param grouping \code{[Deprecated]} Use \code{by} instead.
#'
#' @param labelcex Numeric. Manual control for font size of individual labels. 
#' If \code{NULL} (default), uses a dynamic font size that adjusts automatically 
#' based on matrix dimensions (range 0.2–0.7). Labels are hidden automatically
#' when N > VISMAT_LABEL_MAX (default 500).
#'
#' @param ... Additional arguments passed to the underlying plotting function:
#' \itemize{
#'   \item Heatmap uses \code{\link[lattice]{levelplot}}: can set \code{main},
#'         \code{xlab}, \code{ylab}, \code{col.regions}, \code{colorkey},
#'         \code{scales}, etc.
#'   \item Histogram uses \code{\link[lattice]{histogram}}: can set \code{main},
#'         \code{xlab}, \code{ylab}, \code{nint} (number of bins), etc.
#' }
#'
#' @details
#' \subsection{Compact Matrix Handling}{
#' When \code{mat} is a compact \code{pedmat} object (created with
#' \code{pedmat(..., compact = TRUE)}):
#' \itemize{
#'   \item \strong{With \code{by}}: Group-level mean relationships are computed
#'         algebraically from the K×K compact matrix, including a correction for
#'         sibling off-diagonal values. This avoids expanding to the full N×N
#'         matrix, making family-level or generation-level visualization feasible
#'         even for pedigrees with hundreds of thousands of individuals.
#'   \item \strong{Without \code{by}, N > VISMAT_EXPAND_MAX (5 000)}: The compact K\eqn{\times}K matrix
#'         is plotted directly using representative individuals. Labels show the
#'         number of individuals each representative stands for, e.g.,
#'         \code{"ID (\u00d7350)"}. This avoids memory-intensive full expansion.
#'   \item \strong{Without \code{by}, N \eqn{\le} 5 000}: The compact matrix is
#'         expanded via \code{\link{expand_pedmat}} to restore full dimensions.
#' }
#' }
#'
#' \subsection{Heatmap}{
#' \itemize{
#'   \item Uses a Nature Genetics style color palette (white to orange to red to dark red).
#'   \item Hierarchical clustering reordering (Ward.D2) is enabled by default.
#'   \item Grid lines shown when N \eqn{\le} VISMAT_GRID_MAX (100);
#'         labels shown when N \eqn{\le} VISMAT_LABEL_MAX (500).
#'   \item \code{mat[1,1]} is displayed at the top-left corner.
#' }
#' }
#'
#' \subsection{Histogram}{
#' \itemize{
#'   \item Shows the distribution of lower-triangular elements (pairwise kinship).
#'   \item X-axis: relationship coefficient values; Y-axis: frequency percentage.
#' }
#' }
#'
#' \subsection{Performance}{
#' The following automatic thresholds are defined as package-internal
#' constants (\code{VISMAT_*}) at the top of \code{R/vismat.R}:
#' \itemize{
#'   \item \code{VISMAT_EXPAND_MAX} (5 000): compact matrices with original
#'         N above this are shown in representative view instead of expanding.
#'   \item \code{VISMAT_REORDER_MAX} (2 000): hierarchical clustering is
#'         automatically skipped.
#'   \item \code{VISMAT_LABEL_MAX} (500): individual labels are hidden.
#'   \item \code{VISMAT_GRID_MAX} (100): cell grid lines are hidden.
#'   \item \code{by} grouping uses vectorized \code{rowsum()} algebra — suitable
#'         for large matrices.
#' }
#' }
#'
#' \subsection{Interpreting Relationship Coefficients}{
#' For the additive relationship matrix A:
#' \itemize{
#'   \item Diagonal elements = 1 + F (F = inbreeding coefficient).
#'   \item Off-diagonal elements = 2 × kinship coefficient.
#'   \item 0: unrelated; 0.25: half-sibs / grandparent–grandchild;
#'         0.5: full-sibs / parent–offspring; 1.0: same individual.
#' }
#' }
#'
#' @return Invisibly returns the \code{lattice} plot object. The plot is 
#' rendered on the current graphics device.
#'
#' @seealso 
#' \code{\link{pedmat}} for computing relationship matrices,
#' \code{\link{expand_pedmat}} for manually restoring compact matrix dimensions,
#' \code{\link{query_relationship}} for querying individual pairs,
#' \code{\link{tidyped}} for tidying pedigree data,
#' \code{\link{visped}} for visualizing pedigree structure graphs,
#' \code{\link[lattice]{levelplot}}, \code{\link[lattice]{histogram}}
#'
#' @examples
#' library(visPedigree)
#' data(small_ped)
#' ped <- tidyped(small_ped)
#'
#' # ============================================================
#' # Basic Usage
#' # ============================================================
#'
#' # Method 1: from tidyped object (auto-computes A)
#' vismat(ped)
#'
#' # Method 2: from pedmat object
#' A <- pedmat(ped)
#' vismat(A)
#'
#' # Method 3: from plain matrix
#' vismat(as.matrix(A))
#'
#' # ============================================================
#' # Compact Pedigree (auto-expanded before plotting)
#' # ============================================================
#'
#' # For pedigrees with large full-sib families, compute a compact matrix
#' # first for efficiency, then pass directly to vismat() — it automatically
#' # expands back to full dimensions.
#' A_compact <- pedmat(ped, compact = TRUE)
#' vismat(A_compact)   # prints: "Expanding compact matrix (N -> M individuals)"
#'
#' # For very large pedigrees, aggregate to a group-level view instead
#' vismat(A, ped = ped, by = "Gen",
#'        main = "Mean Relationship Between Generations")
#'
#' # ============================================================
#' # Heatmap Customization
#' # ============================================================
#'
#' # Custom title and axis labels
#' vismat(A, main = "Additive Relationship Matrix",
#'        xlab = "Individual", ylab = "Individual")
#'
#' # Preserve original pedigree order (no clustering)
#' vismat(A, reorder = FALSE)
#'
#' # Custom label font size
#' vismat(A, labelcex = 0.5)
#'
#' # Custom color palette (blue-white-red)
#' vismat(A, col.regions = colorRampPalette(c("blue", "white", "red"))(100))
#'
#' # ============================================================
#' # Display a Subset of Individuals
#' # ============================================================
#'
#' target_ids <- rownames(A)[1:8]
#' vismat(A, ids = target_ids)
#'
#' # ============================================================
#' # Histogram of Relationship Coefficients
#' # ============================================================
#'
#' vismat(A, type = "histogram")
#' vismat(A, type = "histogram", nint = 30)
#'
#' # ============================================================
#' # Group-level Aggregation
#' # ============================================================
#'
#' # Group by generation
#' vismat(A, ped = ped, by = "Gen",
#'        main = "Mean Relationship Between Generations")
#'
#' # Group by full-sib family (founders without a family are excluded)
#' vismat(A, ped = ped, by = "Family")
#'
#' # ============================================================
#' # Other Relationship Matrices
#' # ============================================================
#'
#' # Dominance relationship matrix
#' D <- pedmat(ped, method = "D")
#' vismat(D, main = "Dominance Relationship Matrix")
#' 
#' @export
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom stats as.dist hclust
#' @importFrom lattice levelplot panel.levelplot panel.abline histogram
#' @importFrom data.table as.data.table
vismat <- function(mat, ped = NULL, type = "heatmap", ids = NULL, reorder = TRUE, by = NULL, grouping = NULL, labelcex = NULL, ...) {
  # Backward-compatible deprecation for 'grouping'
  if (!is.null(grouping)) {
    warning("'grouping' is deprecated in vismat(), use 'by' instead.", call. = FALSE)
    if (is.null(by)) by <- grouping
  }
  
  # Track original by label for axis labels (may be consumed early in compact path)
  by_label <- by
  grouping_main <- NULL
  
  # 0a. Extract ped from pedmat object if available
  is_pedmat <- inherits(mat, "pedmat") || !is.null(attr(mat, "pedmat_S4"))
  if (is_pedmat) {
    # Check for inverse matrices - not supported for visualization
    ci <- attr(mat, "call_info")
    if (!is.null(ci) && ci$method[1] %in% c("Ainv", "Dinv", "AAinv")) {
      stop(sprintf(
        paste0("vismat() does not support inverse matrices (%s). ",
               "Inverse matrix elements do not represent meaningful relationship ",
               "coefficients for visualization."),
        ci$method[1]
      ), call. = FALSE)
    }
    
    # Handle compact matrices
    if (isTRUE(ci$compact)) {
      if (!is.null(by)) {
        # Fast path: aggregate directly from K×K compact matrix → G×G
        if (is.null(ped)) ped <- ci$ped_original
        
        # If ids not specified, default to all individuals
        focal_ids <- if (!is.null(ids)) as.character(ids) else ci$ped_original$Ind
        
        agg_mat <- aggregate_compact_by_group(mat, focal_ids, by, ped)
        
        # Replace mat with the G×G result; disable downstream by/ids processing
        mat <- agg_mat
        grouping_main <- sprintf("Grouped Relationship Heatmap (%s)", by)
        by  <- NULL
        ids <- NULL
        
        # Strip pedmat class
        if (inherits(mat, "pedmat")) {
          class(mat) <- setdiff(class(mat), "pedmat")
        }
      } else {
        # No grouping: decide between expand vs compact representative view
        if (ci$n_original > VISMAT_EXPAND_MAX) {
          # Large pedigree: use compact K×K representative view directly
          message(sprintf(
            paste0("Using compact representative view (%d representatives ",
                   "for %d individuals). Use 'by' for group-level aggregation."),
            ci$n_compact, ci$n_original
          ))
          if (is.null(ped)) ped <- ci$ped_original
          cmap <- attr(mat, "compact_map")
          
          # Build informative labels: "ID (×n)" for representatives with siblings
          rep_rows <- cmap[cmap$IsRepresentative == TRUE, ]
          rep_rows <- rep_rows[match(rownames(mat), rep_rows$Ind), ]
          new_labels <- ifelse(
            rep_rows$FamilySize > 1,
            sprintf("%s (\u00d7%d)", rep_rows$Ind, rep_rows$FamilySize),
            rep_rows$Ind
          )
          
          # Strip pedmat class and apply labels
          if (inherits(mat, "pedmat")) {
            class(mat) <- setdiff(class(mat), "pedmat")
          }
          rownames(mat) <- new_labels
          colnames(mat) <- new_labels
          
          # Map ids to representative IDs if provided
          if (!is.null(ids)) {
            ids_char <- as.character(ids)
            rep_map <- cmap[cmap$Ind %chin% ids_char, ]
            rep_ids <- unique(rep_map$RepInd)
            rep_labels <- new_labels[match(rep_ids, rep_rows$Ind)]
            rep_labels <- rep_labels[!is.na(rep_labels)]
            if (length(rep_labels) == 0) {
              stop("None of the specified 'ids' were found in the compact matrix.")
            }
            ids <- rep_labels
          }
        } else {
          # Small enough to expand to full N×N
          message(sprintf(
            "Expanding compact matrix (%d -> %d individuals) for visualization.",
            ci$n_compact, ci$n_original
          ))
          if (is.null(ped)) ped <- ci$ped_original
          mat <- expand_pedmat(mat)
          # Strip pedmat class to get raw matrix
          if (inherits(mat, "pedmat")) {
            class(mat) <- setdiff(class(mat), "pedmat")
          }
        }
      }
    } else {
      if (is.null(ped)) ped <- attr(mat, "ped")
      # Strip pedmat class to get raw matrix
      if (inherits(mat, "pedmat")) {
        class(mat) <- setdiff(class(mat), "pedmat")
      }
    }
  }
  
  # 0b. If input is a tidyped object, calculate A first
  if (is_tidyped(mat)) {
    ped <- mat
    res <- pedmat(ped, method = "A")
    # Now res is a pure matrix with pedmat class
    if (is.null(ped)) {
      ped <- attr(res, "ped")
    }
    class(res) <- setdiff(class(res), "pedmat")
    mat <- res
  }

  # 1. Handle list input from pedmat
  if (is.list(mat) && !is.matrix(mat) && !inherits(mat, "Matrix")) {
    # Reject inverse matrices - primarily via list element names
    # Note: Due to S4 Matrix limitations, unnamed/custom-named lists may bypass this check
    inverse_methods <- c("Ainv", "Dinv", "AAinv")
    
    # Check list element names
    has_inverse_name <- any(names(mat) %in% inverse_methods)
    
    if (has_inverse_name) {
      stop("vismat() does not support inverse matrices (Ainv/Dinv/AAinv). ",
           "Inverse matrix elements do not represent meaningful relationship coefficients. ",
           "Please use pedmat() with method='A', 'D', or 'AA' instead.")
    }
    
    preferred <- c("A", "D", "AA")
    found <- intersect(preferred, names(mat))
    if (length(found) > 0) {
      target_name <- found[1]
      if (length(mat) > 1) {
        message(sprintf("Found multiple matrices in list. Using '%s'.", target_name))
      }
      mat <- mat[[target_name]]
    } else {
      is_mat <- sapply(mat, function(x) is.matrix(x) || inherits(x, "Matrix"))
      if (any(is_mat)) {
        mat <- mat[[which(is_mat)[1]]]
      } else {
        found_types <- sapply(mat, function(x) class(x)[1])
        stop(sprintf("Input list does not contain any valid matrices. Found: %s", paste(found_types, collapse = ", ")))
      }
    }
  }

  # 2. Validation
  if (!is.matrix(mat) && !inherits(mat, "Matrix")) {
    stop("Input 'mat' must be a matrix or a Matrix object.")
  }

  # 3. Subset by IDs if requested
  if (!is.null(ids)) {
    all_names <- rownames(mat)
    if (is.null(all_names)) {
      stop("Matrix has no rownames. Subset by 'ids' is not possible.")
    }
    valid_ids <- intersect(ids, all_names)
    if (length(valid_ids) == 0) {
      stop("None of the specified 'ids' were found in the matrix.")
    }
    mat <- mat[valid_ids, valid_ids, drop = FALSE]
  }

  # 3b. Reorder for clustering (if requested and type is heatmap)
  if (reorder && type == "heatmap") {
    n <- nrow(mat)
    if (n > VISMAT_REORDER_MAX) {
      warning(sprintf("Matrix too large for reordering (N > %d). Skipping clustering.", VISMAT_REORDER_MAX))
      reorder <- FALSE  # Mark as not reordered for downstream logic
    } else {
      # CRITICAL: We use dist() on the matrix rows (Relationship Profile Distance).
      # Siblings share almost identical relationship profiles with the whole population,
      # making their Euclidean distance much smaller than the Parent-Child distance.
      d_mat <- as.matrix(mat)
      dist_obj <- stats::dist(d_mat)
      hc <- stats::hclust(dist_obj, method = "ward.D2")
      mat <- mat[hc$order, hc$order, drop = FALSE]
    }
  }

  # 4. Handle Grouping (Aggregate)
  if (!is.null(by)) {
    if (is.null(ped)) {
      stop("'ped' must be provided when using 'by'.")
    }

    # Ensure ped is a data.table and includes the required column
    ped_dt <- data.table::as.data.table(ped)
    if (!by %in% names(ped_dt)) {
      stop(sprintf("Column '%s' not found in pedigree.", by))
    }

    # Match matrix IDs to groups
    # tidyped objects always use "Ind" as the ID column
    id_col <- "Ind"
    if (!id_col %in% names(ped_dt)) {
      # Fallback to the first column if "Ind" is missing (unlikely for tidyped)
      id_col <- names(ped_dt)[1]
    }
    mapping <- ped_dt[, .(id = as.character(get(id_col)), grp = get(by))]

    # Get group for each matrix row
    mat_ids <- rownames(mat)
    if (is.null(mat_ids)) stop("Matrix must have row names for grouping.")

    # Match matrix IDs to pedigree
    match_idx <- match(mat_ids, mapping$id)
    if (any(is.na(match_idx))) {
      warning("Some individuals in matrix not found in pedigree. They will be assigned to 'Unknown' group.")
    }
    
    # Extract groups (may contain NA if by column has NA values)
    mat_grps <- mapping[match_idx, grp]
    
    # Handle NA groups
    if (any(is.na(mat_grps))) {
      n_na <- sum(is.na(mat_grps))
      
      # For Family grouping, exclude NA individuals (founders without family)
      if (by == "Family") {
        na_idx <- which(is.na(mat_grps))
        na_ids <- mat_ids[na_idx]
        message(sprintf(
          "Note: Excluding %d founder(s) with no family assignment: %s%s",
          n_na, 
          paste(head(na_ids, 5), collapse = ", "),
          if (n_na > 5) sprintf(" (and %d more)", n_na - 5) else ""
        ))
        
        # Remove NA individuals from matrix and by groups
        mat <- mat[-na_idx, -na_idx, drop = FALSE]
        mat_grps <- mat_grps[-na_idx]
        mat_ids <- mat_ids[-na_idx]
      } else {
        # For other by columns, assign to 'Unknown' group
        message(sprintf("Note: %d individual(s) have NA in '%s' column. Assigning to 'Unknown' group.", 
                        n_na, by))
        mat_grps[is.na(mat_grps)] <- "Unknown"
      }
    }

    # Aggregate: Mean relationship between groups
    mat_grps <- as.character(mat_grps)
    grps_unique <- sort(unique(mat_grps))
    n_grp <- length(grps_unique)

    # Check if this is an inverse matrix (averaging inverse elements is usually not useful)
    # Note: Detection based on common naming patterns in row/col names or attributes
    mat_name <- if (!is.null(attr(mat, "method"))) attr(mat, "method") else ""
    is_inv <- grepl("inv", mat_name, ignore.case = TRUE) || 
              grepl("inv", paste(rownames(mat), collapse = ""), ignore.case = TRUE)
    if (is_inv) {
      warning("Aggregating an inverse relationship matrix. Mean values may not be biologically meaningful.")
    }

    message(sprintf("Aggregating %d individuals into %d groups based on '%s'...", nrow(mat), n_grp, by))

    # Optimized aggregation using matrix algebra instead of nested loops
    # This is O(n * n_grp) instead of O(n_grp^2 * avg_group_size^2)
    
    # Create group membership matrix (n x n_grp indicator matrix)
    grp_factor <- factor(mat_grps, levels = grps_unique)
    n <- nrow(mat)
    
    # Count members per group for averaging
    grp_sizes <- tabulate(grp_factor, nbins = n_grp)
    
    # For very large matrices with many groups, use block-wise aggregation
    # to avoid memory issues while being much faster than element-wise loops
    if (n > VISMAT_BLOCK_AGG_MIN_N && n_grp > VISMAT_BLOCK_AGG_MIN_G) {
      message("  Using memory-efficient block aggregation for large matrix...")
      
      # Pre-compute group indices once
      idx_list <- split(seq_len(n), grp_factor)
      
      # Initialize result matrix
      agg_mat <- matrix(0, nrow = n_grp, ncol = n_grp, 
                        dimnames = list(grps_unique, grps_unique))
      
      # Process in row blocks for memory efficiency
      # Aggregate row sums first, then column sums
      for (i in seq_len(n_grp)) {
        idx_i <- idx_list[[i]]
        n_i <- length(idx_i)
        
        # Extract rows for group i (sparse-friendly)
        rows_i <- mat[idx_i, , drop = FALSE]
        
        # Sum across columns for each target group
        for (j in i:n_grp) {
          idx_j <- idx_list[[j]]
          n_j <- length(idx_j)
          
          # Extract the subblock and compute sum (faster than mean on subblocks)
          sub_sum <- sum(rows_i[, idx_j, drop = FALSE])
          agg_mat[i, j] <- sub_sum / (n_i * n_j)
          if (i != j) agg_mat[j, i] <- agg_mat[i, j]
        }
      }
    } else {
      # For smaller matrices, use direct vectorized approach
      # Convert to dense if sparse (for small matrices this is fine)
      mat_dense <- as.matrix(mat)
      
      # Use rowsum for fast aggregation: first sum rows within groups
      # rowsum(mat, group) gives sum of rows for each group
      row_agg <- rowsum(mat_dense, grp_factor, reorder = TRUE)  # n_grp x n
      
      # Then sum columns within groups
      agg_sum <- rowsum(t(row_agg), grp_factor, reorder = TRUE)  # n_grp x n_grp
      agg_sum <- t(agg_sum)
      
      # Compute mean by dividing by group size products
      size_mat <- outer(grp_sizes, grp_sizes)
      agg_mat <- agg_sum / size_mat
      
      # Ensure symmetry and proper naming
      agg_mat <- (agg_mat + t(agg_mat)) / 2
      dimnames(agg_mat) <- list(grps_unique, grps_unique)
    }
    
    # Ensure all values are finite
    agg_mat[!is.finite(agg_mat)] <- 0
    
    mat <- agg_mat
    # Note: User's reorder setting is preserved for grouped matrices

    # Set default titles for later use in heatmap
    grouping_main <- sprintf("Grouped Relationship Heatmap (%s)", by)
  } else {
    grouping_main <- NULL
  }

  # 5. Rendering
  if (type == "heatmap") {
    dots <- list(...)
    if (!"main" %in% names(dots)) {
      dots$main <- if(!is.null(grouping_main)) grouping_main else "Relationship Matrix Heatmap"
    }
    if (!"xlab" %in% names(dots)) dots$xlab <- if(!is.null(by_label)) by_label else (if(reorder) "Clustered Individuals" else "Individuals")
    if (!"ylab" %in% names(dots)) dots$ylab <- if(!is.null(by_label)) by_label else (if(reorder) "Clustered Individuals" else "Individuals")
    
    nature_genetics_palette <- grDevices::colorRampPalette(c(
      "#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84", 
      "#fc8d59", "#ef6548", "#d7301f", "#b30000", "#7f0000"
    ))(100)
    if (!"col.regions" %in% names(dots)) {
      dots$col.regions <- nature_genetics_palette
    }

    # Add a color key (legend) by default
    if (!"colorkey" %in% names(dots)) {
      dots$colorkey <- TRUE
    }

    # IMPORTANT: Matrix::image() from the Matrix package handles dispatch in a way 
    # that often conflicts with lattice arguments for base matrices.
    # We use lattice::levelplot directly to ensure full grid (incl. 0s) and borders.
    mat_to_plot <- as.matrix(mat)
    
    # Final safety check for lattice
    if (!any(is.finite(mat_to_plot))) {
      stop("The relationship matrix contains no finite values to plot.")
    }
    
    # Set up scales for individual labels
    # We want Individual 1 at Top-Left, so Y-axis labels must be reversed.
    # Show labels up to VISMAT_LABEL_MAX with dynamic font scaling.
    if (nrow(mat) <= VISMAT_LABEL_MAX && !"scales" %in% names(dots)) {
      # Use manual labelcex if provided, otherwise dynamic font size
      use_cex <- if(!is.null(labelcex)) labelcex else max(0.2, 0.7 * (1 - (nrow(mat) - 20) / 600))
      
      dots$scales <- list(
        x = list(at = 1:nrow(mat), labels = rownames(mat), rot = 90, cex = use_cex),
        y = list(at = 1:ncol(mat), labels = rev(colnames(mat)), cex = use_cex)
      )
    }

    # Border/Grid logic: Use panel function to draw white grid lines.
    # Limit grid to VISMAT_GRID_MAX to avoid visual clutter.
    if (nrow(mat) <= VISMAT_GRID_MAX && !"panel" %in% names(dots)) {
      dots$panel <- function(x, y, z, ...) {
        lattice::panel.levelplot(x, y, z, ...)
        lattice::panel.abline(h = (0:nrow(mat)) + 0.5, v = (0:nrow(mat)) + 0.5, col = "white", lwd = 0.4)
      }
    }

    # Execute plotting using lattice::levelplot
    # Transpose and reverse rows to make mat[1,1] appear at Top-Left (1, n)
    p <- do.call(lattice::levelplot, c(list(x = t(mat_to_plot[rev(seq_len(nrow(mat_to_plot))), , drop = FALSE])), dots))
    
    print(p)
    return(invisible(p))
    
  } else if (type == "histogram") {
    # Histogram of relationship coefficients
    # Useful for checking the range of inbreeding/kinship
    dots <- list(...)
    
    # Convert to regular matrix if sparse (Matrix package)
    # This ensures lower.tri() works correctly
    if (inherits(mat, "Matrix")) {
      mat <- as.matrix(mat)
    }
    
    vals <- as.numeric(mat[lower.tri(mat)])
    
    if (!"main" %in% names(dots)) dots$main <- "Distribution of Relationship Coefficients (Kinship)"
    if (!"xlab" %in% names(dots)) dots$xlab <- "Relationship Value"
    if (!"ylab" %in% names(dots)) dots$ylab <- "Frequency (%)"
    
    p <- do.call(lattice::histogram, c(list(x = ~vals), dots))
    print(p)
    return(invisible(p))

  } else {
    stop(sprintf("Visualization type '%s' is not supported. Use 'heatmap' or 'histogram'.", type))
  }
}


# --------------------------------------------------------------------------
# Internal: compute G×G group-mean matrix directly from K×K compact matrix
# --------------------------------------------------------------------------
# This avoids expanding the compact matrix to N×N (which may be infeasible
# when N > ~50k) and instead computes the group-level aggregation
# algebraically.
#
# Algorithm
# ---------
# Let A_c be the K×K compact matrix and W an K×G weight matrix where
# W[k, g] = number of individuals in group g whose representative is k.
#
# Step 1  Raw aggregate:   S_raw = W^T %*% A_c %*% W
# Step 2  Sibling correction: for each representative k that stands for
#         n_k > 1 siblings, the off-diagonal sibling value (sib_val_k)
#         differs from the diagonal value A_c[k,k]. The correction is
#         accumulated per group pair.
# Step 3  Mean:            M[g,h] = S_corrected[g,h] / (n_g * n_h)
#
# @param mat       A compact pedmat object (dsCMatrix or similar).
# @param ids       Character vector of focal individual IDs.
# @param by        Character, name of column in ped to group by.
# @param ped       A tidyped data.table with at least Ind and the `by` column.
# @return A dense G×G matrix of group means, with group labels as dimnames.
aggregate_compact_by_group <- function(mat, ids, by, ped) {
  ci   <- attr(mat, "call_info")
  cmap <- data.table::copy(attr(mat, "compact_map"))
  primary_method <- ci$method[1]

  # Strip pedmat class for raw matrix ops
  A_c <- mat
  if (inherits(A_c, "pedmat")) {
    class(A_c) <- setdiff(class(A_c), "pedmat")
  }

  # --- Map focal ids to groups ---
  ped_dt <- data.table::as.data.table(ped)
  focal  <- ped_dt[ped_dt$Ind %chin% ids, .(Ind, grp = as.character(get(by)))]

  # Exclude NA groups (e.g. founders without family)
  na_mask <- is.na(focal$grp)
  if (any(na_mask)) {
    n_na <- sum(na_mask)
    na_ids <- focal$Ind[na_mask]
    if (by == "Family") {
      message(sprintf(
        "Note: Excluding %d founder(s) with no family assignment: %s%s",
        n_na, paste(head(na_ids, 5), collapse = ", "),
        if (n_na > 5) sprintf(" (and %d more)", n_na - 5) else ""
      ))
    } else {
      message(sprintf(
        "Note: %d individual(s) have NA in '%s' column. Excluding from aggregation.",
        n_na, by
      ))
    }
    focal <- focal[!na_mask]
  }

  # Attach representative index from compact_map
  focal <- merge(focal, cmap[, .(Ind, RepIndNum)], by = "Ind", all.x = TRUE)
  if (anyNA(focal$RepIndNum)) {
    stop("Some focal individuals are not present in the compact_map.")
  }

  grps_unique <- sort(unique(focal$grp))
  n_grp <- length(grps_unique)
  K <- nrow(A_c)

  message(sprintf("Aggregating %d individuals into %d groups based on '%s'...",
                  nrow(focal), n_grp, by))

  # --- Build K×G weight matrix W ---
  # W[k, g] = count of focal individuals with RepIndNum == k in group g
  grp_factor <- factor(focal$grp, levels = grps_unique)
  W <- matrix(0L, nrow = K, ncol = n_grp)
  for (i in seq_len(nrow(focal))) {
    W[focal$RepIndNum[i], as.integer(grp_factor[i])] <-
      W[focal$RepIndNum[i], as.integer(grp_factor[i])] + 1L
  }

  # Group sizes
  grp_sizes <- colSums(W)  # length n_grp

  # --- Step 1: raw aggregate S = W^T A_c W ---
  # For sparse A_c, A_c %*% W is efficient
  AcW <- A_c %*% W          # K × G (sparse × dense → dense)
  S   <- crossprod(W, AcW)   # G × G  (W^T (K×G)^T, but crossprod does W^T %*% AcW)

  # --- Step 2: sibling correction ---
  # For each representative k that maps >1 focal individuals, the compact
  # matrix stores the *self* relationship (diagonal) for both the diagonal
  # AND the off-diagonal of the sibling block. The true off-diagonal is
  # sib_val_k, computed from parent relationships.
  # Correction delta per group pair:
  #   g != h:  Δ_gh += W[k,g] * W[k,h] * (sib_val_k - A_c[k,k])
  #   g == h:  Δ_gg += W[k,g] * (W[k,g] - 1) * (sib_val_k - A_c[k,k])
  if (primary_method %in% c("A", "D", "AA")) {
    # A matrix for parent lookups
    A_mat <- if (primary_method == "A") A_c else attr(mat, "A_intermediate")

    if (!is.null(A_mat)) {
      # Find representatives with >1 focal individual mapped
      row_totals <- rowSums(W)
      reps_with_mult <- which(row_totals > 1)

      if (length(reps_with_mult) > 0) {
        for (k in reps_with_mult) {
          s_idx <- cmap$SireNum[cmap$RepIndNum == k][1]
          d_idx <- cmap$DamNum[cmap$RepIndNum == k][1]
          sib_val  <- calc_sib_offdiag(A_mat, s_idx, d_idx, primary_method)
          diag_val <- A_c[k, k]
          delta    <- sib_val - diag_val

          if (abs(delta) < .Machine$double.eps) next

          w_k <- W[k, ]  # length n_grp
          # Off-diagonal correction: w_k[g] * w_k[h] * delta for g != h
          # Diagonal correction:     w_k[g] * (w_k[g] - 1) * delta
          # Combined: outer(w_k, w_k) * delta, then subtract w_k * delta on diag
          # (because diagonal of outer = w_k^2, but we want w_k*(w_k-1))
          correction <- outer(w_k, w_k) * delta
          diag(correction) <- w_k * (w_k - 1) * delta
          S <- S + correction
        }
      }
    }
  }

  # --- Step 3: mean ---
  size_mat <- outer(grp_sizes, grp_sizes)
  agg_mat  <- as.matrix(S) / size_mat
  agg_mat  <- (agg_mat + t(agg_mat)) / 2  # enforce symmetry
  agg_mat[!is.finite(agg_mat)] <- 0
  dimnames(agg_mat) <- list(grps_unique, grps_unique)

  agg_mat
}
