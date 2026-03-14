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
#'   \item A \code{pedmat} object returned by \code{\link{pedmat}}
#'   \item A named list containing matrices (preferring A, D, AA)
#'   \item A \code{\link{tidyped}} object (automatically calculates additive relationship matrix A)
#'   \item A standard \code{matrix} or \code{Matrix} object
#' }
#' \strong{Note}: Inverse matrices (Ainv, Dinv, AAinv) are not supported for 
#' visualization because their elements do not represent meaningful relationship
#' coefficients.
#'
#' @param ped Optional. A tidied pedigree object (\code{tidyped}), used for 
#' extracting labels or grouping information. Required when using the 
#' \code{by} parameter. If \code{mat} is a \code{pedmat} object, 
#' the pedigree can be automatically extracted from its attributes.
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
#' Automatically skipped for large matrices (N > 2000) to improve performance.
#' 
#' \strong{Clustering principle}: Based on relationship profile distance (Euclidean).
#' Full-sibs have nearly identical relationship profiles with the population,
#' so they cluster tightly together.
#'
#' @param by Optional. Column name in \code{ped} to group by (e.g., 
#' \code{"Family"}, \code{"Gen"}, \code{"Year"}). When grouping is enabled:
#' \itemize{
#'   \item Individual-level matrix is aggregated to group-level matrix 
#'         (computing mean relationship coefficients between groups)
#'   \item For \code{"Family"} grouping, founders without family assignment are excluded
#'   \item For other grouping columns, NA values are assigned to \code{"Unknown"} group
#' }
#' This is useful for analyzing the structure of large populations.
#'
#' @param grouping \code{[Deprecated]} Use \code{by} instead.
#'
#' @param labelcex Numeric. Manual control for font size of individual labels. 
#' If \code{NULL} (default), uses dynamic font size that adjusts automatically 
#' based on matrix dimensions (range 0.2-0.7). For matrices with more than 500 
#' individuals, labels are automatically hidden.
#'
#' @param ... Additional arguments passed to the plotting function:
#' \itemize{
#'   \item Heatmap uses \code{\link[lattice]{levelplot}}: can set \code{main},
#'         \code{xlab}, \code{ylab}, \code{col.regions}, \code{colorkey},
#'         \code{scales}, etc.
#'   \item Histogram uses \code{\link[lattice]{histogram}}: can set \code{main},
#'         \code{xlab}, \code{ylab}, \code{nint} (number of bins), etc.
#' }
#'
#' @details
#' \subsection{Visualization Types}{
#' 
#' \strong{Heatmap}:
#' \itemize{
#'   \item Uses Nature Genetics style color palette (white to orange to red to dark red)
#'   \item Hierarchical clustering reordering is enabled by default to group similar individuals
#'   \item Matrix[1,1] is displayed at top-left corner
#'   \item Grid lines shown when N <= 100
#'   \item Individual labels shown when N <= 500
#' }
#' 
#' \strong{Histogram}:
#' \itemize{
#'   \item Shows distribution of lower triangular elements (excluding diagonal)
#'   \item X-axis: relationship coefficient values; Y-axis: frequency percentage
#'   \item Useful for checking population inbreeding levels and kinship structure
#' }
#' }
#'
#' \subsection{Performance Considerations}{
#' \itemize{
#'   \item N > 2000: Hierarchical clustering reordering is automatically skipped
#'   \item N > 500: Individual labels are automatically hidden
#'   \item N > 100: Grid lines are automatically hidden
#'   \item Grouping functionality uses optimized matrix algebra, suitable for large matrices
#' }
#' }
#'
#' \subsection{Interpreting Relationship Coefficients}{
#' For additive relationship matrix A:
#' \itemize{
#'   \item Diagonal elements = 1 + F (where F is the inbreeding coefficient)
#'   \item Off-diagonal elements = 2 x kinship coefficient
#'   \item Value 0: No relationship (unrelated)
#'   \item Value 0.25: Half-sibs or grandparent-grandchild
#'   \item Value 0.5: Full-sibs or parent-offspring
#'   \item Value 1.0: Same individual
#' }
#' }
#'
#' @return Invisibly returns the \code{lattice} plot object. The plot is 
#' generated on the current graphics device.
#'
#' @seealso 
#' \code{\link{pedmat}} for computing relationship matrices
#' \code{\link{tidyped}} for tidying pedigree data
#' \code{\link{visped}} for visualizing pedigree structure graphs
#' \code{\link[lattice]{levelplot}} underlying plotting function for heatmaps
#' \code{\link[lattice]{histogram}} underlying plotting function for histograms
#'
#' @examples
#' # ============================================================
#' # Basic Usage
#' # ============================================================
#' 
#' # Load example data
#' data(simple_ped)
#' ped <- tidyped(simple_ped)
#' 
#' # Method 1: Plot directly from tidyped object (auto-computes A matrix)
#' vismat(ped)
#' 
#' # Method 2: Plot from pedmat object
#' A <- pedmat(ped)
#' vismat(A)
#' 
#' # Method 3: Plot from plain matrix
#' A_dense <- as.matrix(A)
#' vismat(A_dense)
#' 
#' # ============================================================
#' # Heatmap Customization
#' # ============================================================
#' 
#' # Custom title and axis labels
#' vismat(A, main = "Additive Relationship Matrix", xlab = "Individual", ylab = "Individual")
#' 
#' # Disable clustering reorder (preserve original order)
#' vismat(A, reorder = FALSE)
#' 
#' # Custom label font size
#' vismat(A, labelcex = 0.5)
#' 
#' # Custom color palette (blue-white-red)
#' vismat(A, col.regions = colorRampPalette(c("blue", "white", "red"))(100))
#' 
#' # ============================================================
#' # Select Specific Individuals
#' # ============================================================
#' 
#' # Display only a subset of individuals
#' target_ids <- rownames(A)[1:8]
#' vismat(A, ids = target_ids)
#' 
#' # ============================================================
#' # Histogram Visualization
#' # ============================================================
#' 
#' # Relationship coefficient distribution histogram
#' vismat(A, type = "histogram")
#' 
#' # Custom number of bins
#' vismat(A, type = "histogram", nint = 30)
#' 
#' # ============================================================
#' # Group Aggregation (for large populations)
#' # ============================================================
#' 
#' # Group by generation
#' vismat(A, ped = ped, by = "Gen", 
#'        main = "Mean Relationship Between Generations")
#' 
#' # Group by family (if pedigree has Family column)
#' # vismat(A, ped = ped, by = "Family")
#' 
#' # ============================================================
#' # Different Types of Relationship Matrices
#' # ============================================================
#' 
#' # Dominance relationship matrix
#' D <- pedmat(ped, method = "D")
#' vismat(D, main = "Dominance Relationship Matrix")
#' 
#' # Inbreeding coefficient distribution (diagonal elements - 1)
#' A_mat <- as.matrix(A)
#' f_values <- Matrix::diag(A_mat) - 1
#' hist(f_values, main = "Inbreeding Coefficient Distribution", xlab = "Inbreeding (F)")
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
    
    if (is.null(ped)) {
      ped <- attr(mat, "ped")
    }
    # Strip pedmat class to get raw matrix
    if (inherits(mat, "pedmat")) {
      class(mat) <- setdiff(class(mat), "pedmat")
    }
  }
  
  # 0b. If input is a tidyped object, calculate A first
  if (inherits(mat, "tidyped")) {
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
    if (n > 2000) {
      warning("Matrix too large for reordering (N > 2000). Skipping clustering.")
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
    if (n > 10000 && n_grp > 100) {
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
    if (!"xlab" %in% names(dots)) dots$xlab <- if(!is.null(by)) by else (if(reorder) "Clustered Individuals" else "Individuals")
    if (!"ylab" %in% names(dots)) dots$ylab <- if(!is.null(by)) by else (if(reorder) "Clustered Individuals" else "Individuals")
    
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
    # We now allow up to 500 labels with dynamic font scaling.
    if (nrow(mat) <= 500 && !"scales" %in% names(dots)) {
      # Use manual labelcex if provided, otherwise dynamic font size
      use_cex <- if(!is.null(labelcex)) labelcex else max(0.2, 0.7 * (1 - (nrow(mat) - 20) / 600))
      
      dots$scales <- list(
        x = list(at = 1:nrow(mat), labels = rownames(mat), rot = 90, cex = use_cex),
        y = list(at = 1:ncol(mat), labels = rev(colnames(mat)), cex = use_cex)
      )
    }

    # Border/Grid logic: Use panel function to draw white grid lines.
    # Limit grid to N <= 100 to avoid visual clutter.
    if (nrow(mat) <= 100 && !"panel" %in% names(dots)) {
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
