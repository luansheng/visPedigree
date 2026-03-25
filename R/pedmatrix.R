# --- Internal: dgeMatrix fast-coercion cache ---
# Wrap a base matrix as dgeMatrix without any N² symmetry/sparsity scan.
# Used by the sparse=TRUE path of pedmat() when dense storage is appropriate.
# inherits(result, "Matrix") == TRUE; as.matrix() round-trips correctly.
# Class definition is cached on first call to avoid repeated namespace lookups.
.dgeMatrix_cache <- new.env(parent = emptyenv())
.to_dgeMatrix <- function(mat) {
  if (is.null(.dgeMatrix_cache$cls))
    .dgeMatrix_cache$cls <- methods::getClass("dgeMatrix", where = asNamespace("Matrix"))
  dn <- dimnames(mat)
  if (is.null(dn)) dn <- list(NULL, NULL)
  methods::new(.dgeMatrix_cache$cls, x = as.double(mat),
               Dim = as.integer(dim(mat)), Dimnames = dn)
}

#' Genetic Relationship Matrices and Inbreeding Coefficients
#'
#' @description
#' Optimized calculation of additive (A), dominance (D), epistatic (AA) 
#' relationship matrices, their inverses, and inbreeding coefficients (f).
#' Uses Rcpp with Meuwissen & Luo (1992) algorithm for efficient computation.
#'
#' @param ped A tidied pedigree from \code{\link{tidyped}}. Must be a single
#'   pedigree, not a splitped object. For splitped results, use 
#'   \code{pedmat(ped_split$GP1, ...)} to process individual groups.
#' @param method Character, one of:
#' \itemize{
#'   \item \code{"A"}: Additive (numerator) relationship matrix (default)
#'   \item \code{"f"}: Inbreeding coefficients (returns named vector)
#'   \item \code{"Ainv"}: Inverse of A using Henderson's rules (O(n) complexity)
#'   \item \code{"D"}: Dominance relationship matrix
#'   \item \code{"Dinv"}: Inverse of D (requires matrix inversion)
#'   \item \code{"AA"}: Additive-by-additive epistatic matrix (A # A)
#'   \item \code{"AAinv"}: Inverse of AA
#' }
#' @param sparse Logical, if \code{TRUE} returns sparse Matrix (recommended for
#'   large pedigrees). Default is \code{TRUE}.
#' @param invert_method Character, method for matrix inversion (Dinv/AAinv only):
#' \itemize{
#'   \item \code{"auto"}: Auto-detect and use optimal method (default)
#'   \item \code{"sympd"}: Force Cholesky decomposition (faster for SPD matrices)
#'   \item \code{"general"}: Force general LU decomposition
#' }
#' @param threads Integer. Number of OpenMP threads to use. Use 0 to keep
#'   the system/default setting. Currently, multi-threading is explicitly implemented for:
#' \itemize{
#'   \item \code{"D"}: Dominance relationship matrix (significant speedup).
#'   \item \code{"Ainv"}: Inverse of A (only for large pedigrees, n >= 5000).
#' }
#' For \code{"Dinv"}, \code{"AA"}, and \code{"AAinv"}, parallelism depends on the 
#' linked BLAS/LAPACK library (e.g., OpenBLAS, MKL, Accelerate) and is not 
#' controlled by this parameter. Methods \code{"A"} and \code{"f"} are single-threaded.
#' @param compact Logical, if \code{TRUE} compacts full-sibling families by
#'   selecting one representative per family. This dramatically reduces matrix
#'   dimensions for pedigrees with large full-sib groups. See Details.
#' 
#' @details 
#' \strong{API Design:}
#' 
#' Only a single method may be requested per call. This design prevents
#' accidental heavy computations. If multiple matrices are needed, call
#' \code{pedmat()} separately for each method.
#'
#' \strong{Compact Mode (\code{compact = TRUE}):}
#' 
#' Full-siblings share identical relationships with all other individuals.
#' Compact mode exploits this by selecting one representative per full-sib
#' family, dramatically reducing matrix size. For example, a pedigree with
#' 170,000 individuals might compact to 1,800 unique relationship patterns.
#' 
#' Key features:
#' \itemize{
#'   \item \code{\link{query_relationship}(x, id1, id2)}: Query any individual
#'         pair, including merged siblings (automatic lookup)
#'   \item \code{\link{expand_pedmat}(x)}: Restore full matrix dimensions
#'   \item \code{\link{vismat}(x)}: Visualize directly (auto-expands compact)
#' }
#' 
#' \strong{Performance Notes:}
#' \itemize{
#'   \item \strong{Ainv}: O(n) complexity using Henderson's rules. Fast for any size.
#'   \item \strong{Dinv/AAinv}: O(n³) matrix inversion. Practical limits:
#'     \itemize{
#'       \item n < 500: ~10-20 ms
#'       \item n = 1,000: ~40-60 ms  
#'       \item n = 2,000: ~130-150 ms
#'       \item n > 2,000: Consider using \code{compact = TRUE}
#'     }
#'   \item \strong{Memory}: Sparse matrices use ~O(nnz) memory; dense use O(n²)
#' }
#'
#' @return Returns a matrix or vector with S3 class \code{"pedmat"}.
#' 
#' \strong{Object type by method:}
#' \itemize{
#'   \item \code{method="f"}: Named numeric vector of inbreeding coefficients
#'   \item All other methods: Sparse or dense matrix (depending on \code{sparse})
#' }
#' 
#' \strong{S3 Methods:}
#' \itemize{
#'   \item \code{print(x)}: Display matrix with metadata header
#'   \item \code{\link{summary_pedmat}(x)}: Detailed statistics (size, compression, mean, density)
#'   \item \code{dim(x)}, \code{length(x)}, \code{Matrix::diag(x)}, \code{t(x)}: Standard operations
#'   \item \code{x[i, j]}: Subsetting (behaves like underlying matrix)
#'   \item \code{as.matrix(x)}: Convert to base matrix
#' }
#' 
#' \strong{Accessing Metadata (use \code{attr()}, not \code{$}):}
#' \itemize{
#'   \item \code{attr(x, "ped")}: The pedigree used (or compact pedigree if \code{compact=TRUE})
#'   \item \code{attr(x, "ped_compact")}: Compact pedigree (when \code{compact=TRUE})
#'   \item \code{attr(x, "method")}: Calculation method used
#'   \item \code{attr(x, "call_info")}: Full calculation metadata (timing, sizes, etc.)
#'   \item \code{names(attributes(x))}: List all available attributes
#' }
#' 
#' \strong{Additional attributes when \code{compact = TRUE}:}
#' \itemize{
#'   \item \code{attr(x, "compact_map")}: data.table mapping individuals to representatives
#'   \item \code{attr(x, "family_summary")}: data.table summarizing merged families
#'   \item \code{attr(x, "compact_stats")}: Compression statistics (ratio, n_families, etc.)
#' }
#' 
#' @seealso 
#' \code{\link{tidyped}} for preparing pedigree data,
#' \code{\link{query_relationship}} for querying individual pairs,
#' \code{\link{expand_pedmat}} for restoring full dimensions,
#' \code{\link{vismat}} for visualization,
#' \code{\link{inbreed}} for simple inbreeding calculation
#' 
#' @examples
#' # Basic usage with small pedigree
#' library(visPedigree)
#' tped <- tidyped(small_ped)
#' 
#' # --- Additive Relationship Matrix (default) ---
#' A <- pedmat(tped)
#' A["A", "B"]      # Relationship between A and B
#' Matrix::diag(A)  # Diagonal = 1 + F (inbreeding)
#' 
#' # --- Inbreeding Coefficients ---
#' f <- pedmat(tped, method = "f")
#' f["Z1"]  # Inbreeding of individual Z1
#' 
#' # --- Using summary_pedmat() ---
#' summary_pedmat(A)   # Detailed matrix statistics
#' 
#' # --- Accessing Metadata ---
#' attr(A, "ped")              # Original pedigree
#' attr(A, "method")           # "A"
#' names(attributes(A))        # All available attributes
#' 
#' # --- Compact Mode (for large full-sib families) ---
#' A_compact <- pedmat(tped, method = "A", compact = TRUE)
#' 
#' # Query relationships (works for any individual, including merged sibs)
#' query_relationship(A_compact, "Z1", "Z2")  # Full-sibs Z1 and Z2
#' 
#' # View compression statistics
#' attr(A_compact, "compact_stats")
#' attr(A_compact, "family_summary")
#' 
#' # Expand back to full size
#' A_full <- expand_pedmat(A_compact)
#' dim(A_full)  # Original dimensions restored
#' 
#' # --- Inverse Matrices ---
#' Ainv <- pedmat(tped, method = "Ainv")  # Henderson's rules (fast)
#' 
#' # --- Dominance and Epistatic ---
#' D <- pedmat(tped, method = "D")
#' AA <- pedmat(tped, method = "AA")
#' 
#' # --- Visualization (requires display device) ---
#' \dontrun{
#' vismat(A)                       # Heatmap of relationship matrix
#' vismat(A_compact)               # Works with compact matrices
#' vismat(A, by = "Gen")     # Group by generation
#' }
#' 
#' @references
#' Meuwissen, T. H. E., & Luo, Z. (1992). Computing inbreeding coefficients 
#' in large populations. Genetics Selection Evolution, 24(4), 305-313.
#' 
#' Henderson, C. R. (1976). A simple method for computing the inverse of a
#' numerator relationship matrix used in prediction of breeding values.
#' Biometrics, 32(1), 69-83.
#' 
#' @export
pedmat <- function(ped, method = "A", sparse = TRUE, invert_method = "auto", 
                     threads = 0, compact = FALSE) {
  # Check for splitped input - not supported

  if (inherits(ped, "splitped")) {
    stop(paste0(
      "pedmat() does not support 'splitped' objects directly.\n",
      "Use lapply() to process each group separately:\n",
      "  lapply(ped_split[grep('^GP', names(ped_split))], pedmat, method = '", 
      method, "')\n",
      "Or process a single group: pedmat(ped_split$GP1, method = '", method, "')"
    ), call. = FALSE)
  }
  
  # Coerce method to character
  method <- as.character(method)
  
  # Only single method allowed
  if (length(method) != 1) {
    stop(sprintf(
      paste0("Only a single method may be requested per call. ",
             "You requested %d methods: %s. ",
             "Call pedmat() separately for each method."),
      length(method), paste(method, collapse = ", ")
    ), call. = FALSE)
  }
  
  # Validate method
  all_methods <- c("A", "f", "Ainv", "D", "Dinv", "AA", "AAinv")
  if (!method %in% all_methods) {
    stop(sprintf("Invalid method: '%s'. Valid methods: %s", 
                 method, paste(all_methods, collapse = ", ")), call. = FALSE)
  }
  
  # Validate invert_method
  invert_method <- match.arg(invert_method, c("auto", "sympd", "general"))

  # Validate threads
  if (!is.numeric(threads) || length(threads) != 1 || is.na(threads) || threads < 0) {
    stop("'threads' must be a single non-negative integer.")
  }
  threads <- as.integer(threads)
  if (threads > 0) {
    if (!cpp_openmp_available()) {
      warning("OpenMP is not available; 'threads' will be ignored.", call. = FALSE)
    }
    prev_threads <- cpp_set_num_threads(threads)
    on.exit(cpp_set_num_threads(prev_threads), add = TRUE)
  }
  
  # Store original pedigree N for call_info
  n_original <- if (is_tidyped(ped)) nrow(ped) else {
    # Check if it's a data.frame/table before calling nrow
    if (is.data.frame(ped)) nrow(ped) else NA_integer_
  }

  # Handle compact mode
  if (compact) {
    # Store original pedigree for metadata
    ped_original <- if (is_tidyped(ped)) ped else tidyped(ped, addnum = TRUE)
    n_original <- nrow(ped_original)

    # Compact the pedigree
    compact_result <- compact_ped_for_matrix(ped_original)
    
    # For D/AA queries in compact mode, we also need A matrix for full-sibling calculations
    needs_A <- method %in% c("D", "AA")
    
    # Calculate the requested method
    result <- pedmat(
      ped = compact_result$ped_compact,
      method = method,
      sparse = sparse,
      invert_method = invert_method,
      threads = threads,
      compact = FALSE
    )
    
    # If D or AA, extract A matrix from the calculation process (if cached)
    # Note: D and AA calculations internally compute A first, so we can reuse it
    A_matrix <- NULL
    if (needs_A) {
      # Try to extract A from the calculation process (if cached)
      # Otherwise calculate it separately (fallback for safety)
      A_matrix <- tryCatch({
        attr(result, "A_intermediate")
      }, error = function(e) NULL)
      
      if (is.null(A_matrix)) {
        A_matrix <- pedmat(
          ped = compact_result$ped_compact,
          method = "A",
          sparse = sparse,
          compact = FALSE
        )
      }
      # Strip pedmat class for storage
      if (inherits(A_matrix, "pedmat")) {
        class(A_matrix) <- setdiff(class(A_matrix), "pedmat")
      }
    }
    
    # Store A matrix if calculated
    if (!is.null(A_matrix)) {
      attr(result, "A_matrix") <- A_matrix
    }
    
    # Update attributes for compaction
    ci <- attr(result, "call_info")
    ci$compact <- TRUE
    ci$n_original <- n_original
    ci$n_compact <- nrow(compact_result$ped_compact)
    ci$ped_original <- ped_original  # Save for expand_pedmat
    attr(result, "call_info") <- ci
    
    # Add compaction-specific attributes
    attr(result, "compact_map") <- compact_result$compact_map
    attr(result, "family_summary") <- compact_result$family_summary
    attr(result, "compact_stats") <- compact_result$compact_stats
    
    # Mark result as pedmat for compact mode
  if (isS4(result)) {
    attr(result, "pedmat_S4") <- TRUE
  } else {
    class(result) <- c("pedmat", class(result))
  }
  
  return(result)
  }

  # 1. Standardize pedigree
  if (is_tidyped(ped) || all(c("Ind", "Sire", "Dam", "Sex", "Gen", "IndNum", "SireNum", "DamNum") %in% names(ped))) {
    ped <- ensure_complete_tidyped(ped, sprintf("pedmat(method = \"%s\")", method))
  } else if (!is_tidyped(ped)) {
    ped <- tidyped(ped, addnum = TRUE)
  }
  
  # Ensure sorted by IndNum to match matrix indexing
  ped <- data.table::copy(ped)
  data.table::setorder(ped, IndNum)

  sire <- ped$SireNum
  dam <- ped$DamNum
  n <- nrow(ped)
  n_original <- n
  
  # Use a named list to store results
  output <- list()
  
  # Lazy-loading objects to avoid redundant calculations
  # Note: A_dense, D_dense, etc. will fail for 1M individuals. 
  # We only build them if specifically requested in 'method'.
  cache <- list(f_res = NULL, A_dense = NULL, D_dense = NULL, AA_dense = NULL)
  
  get_f_res <- function() {
    if (is.null(cache$f_res)) {
      # Reverting to highly optimized sequential version for f
      # It is very fast and has no parallel overhead or sorting overhead
      cache$f_res <<- cpp_calculate_inbreeding(sire, dam)
    }
    cache$f_res
  }
  
  get_A_dense <- function() {
    if (n > 25000) stop("Pedigree too large (>25k) for dense A matrix. Use 'Ainv' or sparse methods.")
    if (is.null(cache$A_dense)) cache$A_dense <<- cpp_calculate_A(sire, dam)
    cache$A_dense
  }

  get_D_dense <- function() {
    if (n > 25000) stop("Pedigree too large (>25k) for dense D matrix.")
    if (is.null(cache$D_dense)) {
      A_mat <- get_A_dense()
      cache$D_dense <<- cpp_calculate_D(sire, dam, A_mat)
    }
    cache$D_dense
  }

  get_AA_dense <- function() {
    if (n > 25000) stop("Pedigree too large (>25k) for dense AA matrix.")
    if (is.null(cache$AA_dense)) {
      A_mat <- get_A_dense()
      cache$AA_dense <<- cpp_calculate_AA(A_mat)
    }
    cache$AA_dense
  }

  # Process methods (ordered by dependency)
  all_methods <- c("A", "f", "Ainv", "D", "Dinv", "AA", "AAinv")
  active_methods <- intersect(all_methods, method)

  # Check for invalid methods
  invalid <- setdiff(method, all_methods)
  if (length(invalid) > 0) {
    stop(sprintf("Invalid method(s): %s", paste(invalid, collapse = ", ")))
  }

  # Process the single requested method (previously supported multiple)
  for (m in active_methods) {
    if (m == "f") {
      # Prefer getting f from existing A diag, else calculated
      if (!is.null(cache$A_dense)) {
        output$f <- diag(cache$A_dense) - 1.0
      } else {
        output$f <- get_f_res()$f
      }
      
    } else if (m == "A") {
      A_mat <- get_A_dense()
      # A matrices are nearly fully non-zero; converting to dgCMatrix via
      # Matrix::Matrix(dense, sparse=TRUE) scans all N² elements to detect
      # zeros — pure overhead with no storage benefit.  Return a dgeMatrix
      # (dense symmetric) which satisfies inherits(., "Matrix") and coerces
      # to a base matrix via as.matrix() without the N² scan.
      output$A <- if (sparse) .to_dgeMatrix(A_mat) else A_mat
      
    } else if (m == "D") {
      D_mat <- get_D_dense()
      output$D <- if (sparse) .to_dgeMatrix(D_mat) else D_mat
      # Store A matrix for compact mode reuse (on output list, not S4 matrix)
      if (!is.null(cache$A_dense)) {
        output$A_intermediate_D <- cache$A_dense
      }
      
    } else if (m == "Dinv") {
      D_mat <- get_D_dense()
      Dinv <- switch(invert_method,
        "auto" = cpp_invert_auto(D_mat),
        "sympd" = cpp_invert_sympd(D_mat),
        "general" = cpp_invert_dense(D_mat)
      )
      output$Dinv <- if (sparse) .to_dgeMatrix(Dinv) else Dinv

    } else if (m == "AA") {
      AA_mat <- get_AA_dense()
      output$AA <- if (sparse) .to_dgeMatrix(AA_mat) else AA_mat
      # Store A matrix for compact mode reuse (on output list, not S4 matrix)
      if (!is.null(cache$A_dense)) {
        output$A_intermediate_AA <- cache$A_dense
      }

    } else if (m == "AAinv") {
      AA_mat <- get_AA_dense()
      AAinv <- switch(invert_method,
        "auto" = cpp_invert_auto(AA_mat),
        "sympd" = cpp_invert_sympd(AA_mat),
        "general" = cpp_invert_dense(AA_mat)
      )
      output$AAinv <- if (sparse) .to_dgeMatrix(AAinv) else AAinv
      
    } else if (m == "Ainv") {
      f_res <- get_f_res()
      # Pass threads if we want, but currently C++ uses max
      triplets <- cpp_build_ainv_triplets(sire, dam, f_res$dii)
      Ainv <- Matrix::sparseMatrix(
        i = triplets$i, j = triplets$j, x = triplets$v,
        dims = c(n, n), symmetric = TRUE
      )
      output$Ainv <- if (sparse) Ainv else as.matrix(Ainv)
    }
  }

  # 3. Add dimension names to matrices
  ind_names <- ped$Ind
  for (m in names(output)) {
    obj <- output[[m]]
    if (is.matrix(obj) || inherits(obj, "Matrix")) {
      # Square matrices: set both dimensions
      rownames(obj) <- ind_names
      colnames(obj) <- ind_names
      output[[m]] <- obj
    } else if (is.numeric(obj) && length(obj) == n) {
      # Vectors (f)
      names(obj) <- ind_names
      output[[m]] <- obj
    }
  }

  # Construct return object: ALWAYS return pure matrix/vector
  # All metadata stored in attributes
  result <- if (length(active_methods) == 1) {
    output[[active_methods]]
  } else {
    # Multiple methods: return named list
    output
  }
  
  # Transfer A_intermediate from output list to result attributes (if exists)
  # This is for compact mode reuse
  if (length(active_methods) == 1) {
    if (active_methods == "D" && !is.null(output$A_intermediate_D)) {
      attr(result, "A_intermediate") <- output$A_intermediate_D
    } else if (active_methods == "AA" && !is.null(output$A_intermediate_AA)) {
      attr(result, "A_intermediate") <- output$A_intermediate_AA
    }
  }
  
  # Attach call information
  attr(result, "call_info") <- list(
    method = if(length(active_methods) == 1) active_methods else names(output),
    sparse = sparse,
    compact = FALSE,  # This branch is never compact
    n_original = n_original,
    n_compact = n_original,
    timestamp = Sys.time()
  )
  
  # Attach method name for easy access
  attr(result, "method") <- if(length(active_methods) == 1) active_methods else names(output)

  # Attach the pedigree used for calculation
  attr(result, "ped") <- ped
  
  # Store all methods if multiple were calculated
  if (length(active_methods) > 1) {
    attr(result, "mat_all") <- output
  }
  
  # Set S3 class - but only for non-S4 objects
  # For S4 objects (like Matrix), we keep them as-is and add pedmat_S4 marker
  if (isS4(result)) {
    attr(result, "pedmat_S4") <- TRUE
  } else {
    class(result) <- c("pedmat", class(result))
  }
  
  return(result)
}


#' Compact pedigree by merging full siblings for matrix calculation
#' 
#' @description
#' This internal function identifies full siblings (individuals sharing the same
#' sire and dam) and selects one representative per family. This can dramatically
#' reduce memory requirements when calculating relationship matrices for pedigrees
#' with large full-sibling families.
#' 
#' @param ped A tidyped object or pedigree data.
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{ped_compact}: The compacted pedigree (tidyped object)
#'   \item \code{compact_map}: Individual-level mapping table
#'   \item \code{family_summary}: Family-level summary table
#'   \item \code{compact_stats}: Compression statistics
#' }
#' 
#' @keywords internal
compact_ped_for_matrix <- function(ped) {
  # Ensure tidyped format
  if (!is_tidyped(ped)) {
    ped <- tidyped(ped, addnum = TRUE)
  }
  
  ped_dt <- data.table::as.data.table(data.table::copy(ped))
  n_original <- nrow(ped_dt)
  
  # Step 1: Create family identifiers
  # Using numeric keys for uniqueness and efficiency
  ped_dt[, family_key := paste(SireNum, DamNum, sep = "_")]
  ped_dt[SireNum == 0 | DamNum == 0, family_key := NA_character_]
  
  # Step 2: Create human-readable family labels
  ped_dt[!is.na(family_key), family_label := paste0(Sire, "x", Dam)]
  
  # Step 3: Calculate family sizes
  ped_dt[!is.na(family_key), family_size := .N, by = family_key]
  ped_dt[is.na(family_key), family_size := 1L]
  
  # Step 4: Identify individuals who are parents
  # Parents must be kept because their offspring reference them
  all_sires <- unique(ped_dt$Sire[!is.na(ped_dt$Sire) & ped_dt$Sire != ""])
  all_dams <- unique(ped_dt$Dam[!is.na(ped_dt$Dam) & ped_dt$Dam != ""])
  all_parents <- unique(c(all_sires, all_dams))
  
  # Mark who is a parent
  ped_dt[, IsParent := Ind %in% all_parents]
  
  # Step 5: Identify full-sibling families (family_size >= 2)
  fullsib_families <- ped_dt[!is.na(family_key) & family_size >= 2]
  
  if (nrow(fullsib_families) == 0) {
    # No full-sibling families to compact
    ped_compact <- data.table::copy(ped)
    return(list(
      ped_compact = ped_compact,
      compact_map = data.table(
        Ind = ped_dt$Ind,
        IndNum = ped_dt$IndNum,
        RepInd = ped_dt$Ind,
        RepIndNum = ped_dt$IndNum,
        FamilyID = NA_character_,
        SireNum = ped_dt$SireNum,
        DamNum = ped_dt$DamNum,
        Sire = ped_dt$Sire,
        Dam = ped_dt$Dam,
        FamilyLabel = NA_character_,
        FamilySize = 1L,
        IsRepresentative = TRUE,
        IsCompacted = FALSE
      ),
      family_summary = data.table(),
      compact_stats = list(
        n_original = n_original,
        n_compact = n_original,
        n_removed = 0L,
        n_families_compacted = 0L,
        compression_ratio = 1.0,
        memory_saved_pct = 0.0
      )
    ))
  }
  
  # Step 6: For each family, select representative from NON-PARENT members only
  # Parents are NEVER removed and NEVER represent others to avoid identity loss.
  
  compactable_members <- fullsib_families[IsParent == FALSE]
  
  if (nrow(compactable_members) == 0) {
    # No non-parent full-siblings to compact
    ped_compact <- data.table::copy(ped)
    return(list(
      ped_compact = ped_compact,
      compact_map = data.table(
        Ind = ped_dt$Ind,
        IndNum = ped_dt$IndNum,
        RepInd = ped_dt$Ind,
        RepIndNum = ped_dt$IndNum,
        FamilyID = NA_character_,
        SireNum = ped_dt$SireNum,
        DamNum = ped_dt$DamNum,
        Sire = ped_dt$Sire,
        Dam = ped_dt$Dam,
        FamilyLabel = NA_character_,
        FamilySize = 1L,
        IsRepresentative = TRUE,
        IsCompacted = FALSE
      ),
      family_summary = data.table(),
      compact_stats = list(
        n_original = n_original,
        n_compact = n_original,
        n_removed = 0L,
        n_families_compacted = 0L,
        compression_ratio = 1.0,
        memory_saved_pct = 0.0
      )
    ))
  }
  
  representatives <- compactable_members[, {
    list(RepInd = Ind[1],
         RepIndNum = IndNum[1])
  }, by = family_key]
  
  # Step 7: Create family lookup table with FamilyID
  family_lookup <- compactable_members[, .(
    SireNum = first(SireNum),
    DamNum = first(DamNum),
    Sire = first(Sire),
    Dam = first(Dam),
    FamilyLabel = first(family_label),
    FamilySize = .N,
    Gen = first(Gen)
  ), by = family_key]
  
  family_lookup[, FamilyID := sprintf("F%04d", .I)]
  
  # Add representative information
  family_lookup <- representatives[family_lookup, on = "family_key"]
  
  # Step 8: Build compact_map for ALL individuals
  # Map compactable members to their representative, others to themselves
  
  compact_map_fullsibs <- family_lookup[compactable_members, on = "family_key",
    .(Ind = i.Ind,
      OldIndNum = i.IndNum,
      RepInd = RepInd,
      OldRepIndNum = RepIndNum,
      FamilyID = FamilyID,
      SireNum = i.SireNum,
      DamNum = i.DamNum,
      Sire = i.Sire,
      Dam = i.Dam,
      FamilyLabel = FamilyLabel,
      FamilySize = FamilySize,
      IsParent = FALSE,
      IsRepresentative = (i.Ind == RepInd),
      IsCompacted = TRUE)
  ]
  
  # Non-compacted individuals (all parents + unique non-parents)
  non_compacted <- ped_dt[!(Ind %in% compactable_members$Ind)]
  compact_map_others <- data.table(
    Ind = non_compacted$Ind,
    OldIndNum = non_compacted$IndNum,
    RepInd = non_compacted$Ind,
    OldRepIndNum = non_compacted$IndNum,
    FamilyID = NA_character_,
    SireNum = non_compacted$SireNum,
    DamNum = non_compacted$DamNum,
    Sire = non_compacted$Sire,
    Dam = non_compacted$Dam,
    FamilyLabel = NA_character_,
    FamilySize = 1L,
    IsParent = non_compacted$IsParent,
    IsRepresentative = TRUE,
    IsCompacted = FALSE
  )
  
  compact_map <- rbind(compact_map_fullsibs, compact_map_others)
  data.table::setorder(compact_map, OldIndNum)
  
  # Step 9: Create compacted pedigree
  # Remove rule: Non-parent AND non-representative individuals
  # Keep rule: All parents + representatives
  remove_ids <- compact_map[IsParent == FALSE & IsRepresentative == FALSE, OldIndNum]
  ped_compact <- ped_dt[!(IndNum %in% remove_ids)]
  data.table::setorder(ped_compact, IndNum)
  
  # CRITICAL: Renumber IndNum to be consecutive (1, 2, 3, ...)
  # C++ code expects IndNum to be consecutive indices starting from 1
  old_to_new_map <- data.table(
    OldNum = ped_compact$IndNum,
    NewNum = seq_len(nrow(ped_compact))
  )
  
  # Update compact_map: map OldRepIndNum to NewRepIndNum
  compact_map <- old_to_new_map[compact_map, on = c(OldNum = "OldRepIndNum")]
  setnames(compact_map, "NewNum", "RepIndNum")
  compact_map[, OldNum := NULL]
  
  # Update SireNum and DamNum in compact_map to match compact matrix indices
  # First SireNum
  compact_map <- old_to_new_map[compact_map, on = c(OldNum = "SireNum")]
  setnames(compact_map, "NewNum", "SireNum")
  compact_map[is.na(SireNum), SireNum := 0L] # Handle founders / missing from map
  compact_map[, OldNum := NULL]
  
  # Then DamNum
  compact_map <- old_to_new_map[compact_map, on = c(OldNum = "DamNum")]
  setnames(compact_map, "NewNum", "DamNum")
  compact_map[is.na(DamNum), DamNum := 0L]
  compact_map[, OldNum := NULL]
  
  # Also store old IndNum for reference
  setnames(compact_map, "OldIndNum", "IndNum")
  
  # Update ped_compact: renumber IndNum, SireNum, DamNum
  ped_compact[, OldIndNum := IndNum]
  ped_compact <- old_to_new_map[ped_compact, on = c(OldNum = "OldIndNum")]
  ped_compact[, IndNum := NewNum]
  ped_compact[, c("NewNum", "OldNum") := NULL]
  
  # Update SireNum
  old_to_new_sire <- copy(old_to_new_map)
  setnames(old_to_new_sire, c("OldSireNum", "NewSireNum"))
  ped_compact <- old_to_new_sire[ped_compact, on = c(OldSireNum = "SireNum")]
  ped_compact[!is.na(NewSireNum), SireNum := NewSireNum]
  ped_compact[is.na(NewSireNum), SireNum := 0L]
  ped_compact[, c("OldSireNum", "NewSireNum") := NULL]
  
  # Update DamNum
  old_to_new_dam <- copy(old_to_new_map)
  setnames(old_to_new_dam, c("OldDamNum", "NewDamNum"))
  ped_compact <- old_to_new_dam[ped_compact, on = c(OldDamNum = "DamNum")]
  ped_compact[!is.na(NewDamNum), DamNum := NewDamNum]
  ped_compact[is.na(NewDamNum), DamNum := 0L]
  ped_compact[, c("OldDamNum", "NewDamNum") := NULL]
  
  # Remove OldIndNum if it still exists
  if ("OldIndNum" %in% names(ped_compact)) {
    ped_compact[, OldIndNum := NULL]
  }
  
  # Clean up extra columns added during compaction
  extra_cols <- c("family_key", "family_label", "family_size", "IsParent")
  for (col in extra_cols) {
    if (col %in% names(ped_compact)) {
      ped_compact[, (col) := NULL]
    }
  }
  
  # Count before tidyped (for accurate statistics)
  n_before_tidyped <- nrow(ped_compact)
  
  # Calculate family compression counts BEFORE tidyped (while we still have compact_map with removed individuals)
  family_compressed_counts <- compact_map[IsCompacted == TRUE, .(
    NCompressed = sum(IsParent == FALSE & IsRepresentative == FALSE)
  ), by = FamilyID]
  
  # CRITICAL: Use tidyped to ensure ped_compact is complete and properly numbered
  # This will add missing parents, renumber everything correctly, and compute Family field
  ped_compact <- tidyped(ped_compact, addnum = TRUE)
  
  # Verify Family and FamilySize columns exist in ped_compact
  if (!("Family" %in% names(ped_compact))) {
    warning("Family field missing in compact pedigree. This may cause issues with vismat(by = 'Family').")
  }
  
  # After tidyped, rebuild compact_map to include ALL individuals (kept + removed)
  # Save original compact_map information before rebuilding
  original_compact_map <- copy(compact_map)
  
  # Build compact_map for individuals in ped_compact (kept individuals)
  kept_compact_map <- ped_compact[, .(
    Ind,
    IndNum,
    SireNum,
    DamNum,
    Sire,
    Dam
  )]
  
  # Merge with original mapping information
  kept_compact_map <- original_compact_map[kept_compact_map, on = "Ind",
    .(Ind = i.Ind,
      IndNum = i.IndNum,
      SireNum = i.SireNum,
      DamNum = i.DamNum,
      Sire = i.Sire,
      Dam = i.Dam,
      RepInd = RepInd,
      FamilyID = FamilyID,
      FamilyLabel = FamilyLabel,
      FamilySize = FamilySize,
      IsCompacted = IsCompacted,
      IsParent = IsParent,
      IsRepresentative = IsRepresentative)
  ]
  
  # For individuals added by tidyped (not in original), set defaults
  kept_compact_map[is.na(RepInd), `:=`(
    RepInd = Ind,
    FamilyID = NA_character_,
    FamilyLabel = NA_character_,
    FamilySize = 1L,
    IsCompacted = FALSE,
    IsParent = TRUE,  # Added by tidyped means they're parents
    IsRepresentative = TRUE
  )]
  
  # Update RepIndNum based on new numbering in ped_compact
  rep_num_mapping <- ped_compact[, .(RepInd = Ind, RepIndNum = IndNum)]
  kept_compact_map <- rep_num_mapping[kept_compact_map, on = "RepInd"]
  
  # Build compact_map for removed individuals (not in ped_compact)
  removed_inds <- original_compact_map[!(Ind %in% ped_compact$Ind)]
  if (nrow(removed_inds) > 0) {
    # For removed individuals, keep original info but update RepIndNum
    removed_compact_map <- removed_inds[, .(
      Ind,
      IndNum = NA_integer_,  # No IndNum in ped_compact
      RepInd,
      SireNum,
      DamNum,
      Sire,
      Dam,
      FamilyID,
      FamilyLabel,
      FamilySize,
      IsCompacted,
      IsParent,
      IsRepresentative = FALSE  # Removed individuals are never representatives
    )]
    
    # Update RepIndNum for removed individuals
    removed_compact_map <- rep_num_mapping[removed_compact_map, on = "RepInd"]
    
    # Combine kept and removed
    new_compact_map <- rbind(kept_compact_map, removed_compact_map, fill = TRUE)
  } else {
    new_compact_map <- kept_compact_map
  }
  
  # CRITICAL: Re-map SireNum and DamNum to match the FINAL ped_compact numbering
  # because tidyped() might have reordered/added individuals.
  final_num_mapping <- ped_compact[, .(ParentID = Ind, NewParentNum = IndNum)]
  
  # Update SireNum
  new_compact_map[, SireNum := final_num_mapping$NewParentNum[match(Sire, final_num_mapping$ParentID)]]
  new_compact_map[is.na(SireNum), SireNum := 0L]
  
  # Update DamNum
  new_compact_map[, DamNum := final_num_mapping$NewParentNum[match(Dam, final_num_mapping$ParentID)]]
  new_compact_map[is.na(DamNum), DamNum := 0L]
  
  # Sort: kept individuals first (by IndNum), then removed individuals (by Ind)
  new_compact_map[!is.na(IndNum), sort_key := IndNum]
  new_compact_map[is.na(IndNum), sort_key := 1000000 + seq_len(.N)]
  data.table::setorder(new_compact_map, sort_key)
  new_compact_map[, sort_key := NULL]
  
  compact_map <- new_compact_map
  data.table::setorder(compact_map, IndNum)
  
  n_compact <- nrow(ped_compact)
  n_removed <- n_original - n_before_tidyped  # Use pre-tidyped count
  n_families <- nrow(family_lookup)
  
  # Step 9: Calculate statistics
  compression_ratio <- n_compact / n_original
  # Memory saved for dense matrix: (n_original^2 - n_compact^2) / n_original^2
  memory_saved_pct <- (1 - (n_compact / n_original)^2) * 100
  
  # By sex statistics
  sex_stats <- compact_map[, .(
    n_original = .N,
    n_compact = sum(IsRepresentative),
    n_removed = sum(!IsRepresentative)
  ), by = .(Sex = ifelse(is.na(Sire) & is.na(Dam), "founder", 
                         ifelse(Sire != "", "with_parents", "unknown")))]
  
  # Family size distribution
  family_size_breaks <- c(0, 2, 10, 50, 100, Inf)
  family_size_labels <- c("1", "2-10", "11-50", "51-100", "100+")
  family_lookup[, SizeCategory := cut(FamilySize, 
                                        breaks = family_size_breaks,
                                        labels = family_size_labels,
                                        right = TRUE)]
  
  family_size_dist <- family_lookup[, .(
    n_families = .N,
    n_individuals_total = sum(FamilySize),
    n_individuals_removed = sum(FamilySize - 1)
  ), by = SizeCategory]
  
  # Ensure result is not a tidyped object (just a data.table)
  class(family_size_dist) <- c("data.table", "data.frame")
  
  compact_stats <- list(
    n_original = n_original,
    n_compact = n_compact,
    n_removed = n_removed,
    n_families_compacted = n_families,
    compression_ratio = compression_ratio,
    memory_saved_pct = memory_saved_pct,
    by_sex = sex_stats,
    family_size_dist = family_size_dist
  )
  
  # Step 10: Prepare final family_summary (update RepIndNum with new numbering)
  # Match RepInd with the new IndNum from ped_compact
  rep_num_mapping <- ped_compact[, .(RepInd = Ind, NewRepIndNum = IndNum)]
  family_summary <- family_lookup[, .(FamilyID, SireNum, DamNum, Sire, Dam, 
                                      FamilyLabel, FamilySize, RepInd, Gen)]
  family_summary <- rep_num_mapping[family_summary, on = "RepInd"]
  setnames(family_summary, "NewRepIndNum", "RepIndNum")
  

  # Add NCompressed column using pre-calculated counts (before tidyped)
  family_summary <- family_compressed_counts[family_summary, on = "FamilyID"]
  family_summary[is.na(NCompressed), NCompressed := 0L]
  
  # Reorder columns for clarity
  data.table::setcolorder(family_summary, c("FamilyID", "FamilyLabel", "Sire", "Dam",
                                             "SireNum", "DamNum", "FamilySize", 
                                             "NCompressed", "RepInd", "RepIndNum", "Gen"))
  data.table::setorder(family_summary, FamilyID)
  
  # Ensure result is not a tidyped object
  class(family_summary) <- c("data.table", "data.frame")
  class(compact_map) <- c("data.table", "data.frame")
  
  return(list(
    ped_compact = ped_compact,
    compact_map = compact_map,
    family_summary = family_summary,
    compact_stats = compact_stats
  ))
}


#' Query Relationship Coefficients from a Pedigree Matrix
#' 
#' @description
#' Retrieves relationship coefficients between individuals from a pedmat object.
#' For compact matrices, automatically handles lookup of merged full-siblings.
#' 
#' @param x A pedmat object created by \code{\link{pedmat}}.
#' @param id1 Character, first individual ID.
#' @param id2 Character, second individual ID. If \code{NULL}, returns 
#'   the entire row of relationships for \code{id1}.
#' 
#' @return 
#' \itemize{
#'   \item If \code{id2} is provided: numeric value (relationship coefficient)
#'   \item If \code{id2} is \code{NULL}: named numeric vector (id1's row)
#'   \item Returns \code{NA} if individual not found
#' }
#' 
#' @details
#' For compact matrices (\code{compact = TRUE}), this function automatically

#' maps individuals to their family representatives. For methods A, D, and AA,
#' it can compute the correct relationship even between merged full-siblings
#' using the formula:
#' \itemize{
#'   \item \strong{A}: \eqn{a_{ij} = 0.5 * (a_{is} + a_{id})} where s, d are parents
#'   \item \strong{D}: \eqn{d_{ij} = a_{ij}^2} (for full-sibs in same family)
#'   \item \strong{AA}: \eqn{aa_{ij} = a_{ij}^2}
#' }
#' 
#' @note
#' Inverse matrices (Ainv, Dinv, AAinv) are \strong{not supported} because
#' inverse matrix elements do not represent meaningful relationship coefficients.
#' 
#' @seealso \code{\link{pedmat}}, \code{\link{expand_pedmat}}
#' 
#' @examples
#' tped <- tidyped(small_ped)
#' A <- pedmat(tped, method = "A", compact = TRUE)
#' 
#' # Query specific pair
#' query_relationship(A, "A", "B")
#' 
#' # Query merged full-siblings (works with compact)
#' query_relationship(A, "Z1", "Z2")
#' 
#' # Get all relationships for one individual
#' query_relationship(A, "A")
#' 
#' @export
query_relationship <- function(x, id1, id2 = NULL) {
  # Check if x is a pedmat (either S3 or S4 with marker)
  is_pedmat <- inherits(x, "pedmat") || !is.null(attr(x, "pedmat_S4"))
  if (!is_pedmat) {
    stop("x must be a pedmat object")
  }
  
  ci <- attr(x, "call_info")
  primary_method <- ci$method[1]
  
  # Inverse matrices are not supported - their elements don't represent
  # meaningful relationship coefficients
  if (primary_method %in% c("Ainv", "Dinv", "AAinv")) {
    stop(sprintf(
      paste0("query_relationship() does not support inverse matrices (%s). ",
             "Inverse matrix elements do not represent meaningful relationship ",
             "coefficients between individuals."),
      primary_method
    ), call. = FALSE)
  }
  
  # Handle vector case (e.g., inbreeding coefficients "f")
  # For vectors, id2 is ignored - we just look up the value for id1
  if (primary_method == "f" || (is.numeric(x) && is.null(dim(x)))) {
    # Strip class for indexing
    vec <- x
    if (inherits(vec, "pedmat")) {
      class(vec) <- setdiff(class(vec), "pedmat")
    }
    
    if (isTRUE(ci$compact)) {
      compact_map <- attr(x, "compact_map")
      map1 <- compact_map[Ind == id1]
      if (nrow(map1) == 0) return(NA_real_)
      rep1_idx <- map1$RepIndNum[1]
      return(vec[rep1_idx])
    } else {
      # Non-compact: use names directly
      if (id1 %in% names(vec)) {
        return(vec[id1])
      } else {
        return(NA_real_)
      }
    }
  }
  
  # Get raw matrix (already is the raw matrix)
  mat <- x
  if (inherits(mat, "pedmat")) {
    class(mat) <- setdiff(class(mat), "pedmat")
  }
  
  # Check if it's compact
  if (isTRUE(ci$compact)) {
    compact_map <- attr(x, "compact_map")
    
    # Use data.table for fast lookup
    map1 <- compact_map[Ind == id1]
    if (nrow(map1) == 0) return(NA_real_)
    rep1_idx <- map1$RepIndNum[1]
    rep1_ind <- map1$RepInd[1]
    
    if (is.null(id2)) {
      # Return row for representative
      # The returned vector keeps the matrix column names (compact pedigree individuals)
      res <- mat[rep1_idx, ]
      # Names are already set from matrix colnames, no need to modify
      return(res)
    }
    
    map2 <- compact_map[Ind == id2]
    if (nrow(map2) == 0) return(NA_real_)
    rep2_idx <- map2$RepIndNum[1]
    rep2_ind <- map2$RepInd[1]
    
    # The actual individuals are different but they share the same representative
    # (Full-sibling case)
    if (id1 != id2 && rep1_idx == rep2_idx) {
      repS_idx <- map1$SireNum[1]
      repD_idx <- map1$DamNum[1]
      
      # Helper to get relationship (handling missing parents as 0)
      # We need A values even if primary_method is D or AA
      get_A_val <- function(r1, r2) {
        if (r1 <= 0 || r2 <= 0) return(0.0)
        
        # Look for A matrix
        A_mat <- NULL
        if (primary_method == "A") {
            A_mat <- mat
        } else {
            # Check for stored A matrix in attributes
            # Try both "A_matrix" and "A_intermediate" for backward/internal compatibility
            A_mat <- attr(x, "A_intermediate")
            if (is.null(A_mat)) {
              A_mat <- attr(x, "A_matrix")
            }
            
            if (is.null(A_mat)) {
              mat_all <- attr(x, "mat_all")
              if (!is.null(mat_all) && !is.null(mat_all$A)) {
                A_mat <- mat_all$A
              }
            }
        }
        
        if (is.null(A_mat)) return(NA_real_)
        return(A_mat[r1, r2])
      }
      
      A_ss <- get_A_val(repS_idx, repS_idx)
      if (repS_idx <= 0) A_ss <- 1.0 
      A_dd <- get_A_val(repD_idx, repD_idx)
      if (repD_idx <= 0) A_dd <- 1.0
      A_sd <- get_A_val(repS_idx, repD_idx)
      
      
      # Calculate A_ij
      if (repS_idx <= 0 && repD_idx <= 0) {
        res_A <- 0.0
      } else {
        res_A <- 0.25 * (A_ss + A_dd + 2 * A_sd)
      }
      
      if (primary_method == "A") return(res_A)
      
      # For D: 0.25 * (A_SS * A_DD + A_SD^2)
      if (primary_method == "D") return(0.25 * (A_ss * A_dd + A_sd^2))
      
      # For AA: A_ij^2
      if (primary_method == "AA") return(res_A^2)
      
      return(res_A) 
    }
    
    return(mat[rep1_idx, rep2_idx])
  } else {
    # Non-compact: direct indexing
    if (is.null(id2)) return(mat[id1, ])
    return(mat[id1, id2])
  }
}

# --------------------------------------------------------------------------
# Internal helper: sibling off-diagonal value from parent relationships
# --------------------------------------------------------------------------
# Computes the relationship value between full-siblings who share the same
# compact representative. Used by expand_pedmat() and aggregate_compact_by_group().
#
# @param A_mat The compact A (additive) matrix.
# @param sire_idx Integer, row index of the sire in A_mat (0 = unknown).
# @param dam_idx  Integer, row index of the dam in A_mat (0 = unknown).
# @param method   Character, one of "A", "D", "AA".
# @return Numeric scalar: the sibling off-diagonal value.
calc_sib_offdiag <- function(A_mat, sire_idx, dam_idx, method = "A") {
  A_ss <- if (sire_idx > 0) A_mat[sire_idx, sire_idx] else 1.0
  A_dd <- if (dam_idx  > 0) A_mat[dam_idx,  dam_idx]  else 1.0
  A_sd <- if (sire_idx > 0 && dam_idx > 0) A_mat[sire_idx, dam_idx] else 0.0

  res_A_sib <- 0.25 * (A_ss + A_dd + 2 * A_sd)

  switch(method,
    "A"  = res_A_sib,
    "D"  = 0.25 * (A_ss * A_dd + A_sd^2),
    "AA" = res_A_sib^2,
    stop(sprintf("Unsupported method '%s' for sibling off-diagonal.", method))
  )
}

#' Expand a Compact Pedigree Matrix to Full Dimensions
#' 
#' @description
#' Restores a compact pedmat to its original dimensions by mapping each
#' individual to their family representative's values. For non-compact matrices,
#' returns the matrix unchanged.
#' 
#' @param x A pedmat object from \code{\link{pedmat}}.
#' 
#' @return 
#' Matrix or vector with original pedigree dimensions:
#' \itemize{
#'   \item Matrices: Row and column names set to all individual IDs
#'   \item Vectors (e.g., method="f"): Names set to all individual IDs
#' }
#' The result is \strong{not} a pedmat object (S3 class stripped).
#' 
#' @details
#' For compact matrices, full-siblings within the same family will have
#' identical relationship values in the expanded matrix because they shared
#' the same representative during calculation.
#' 
#' @seealso \code{\link{pedmat}}, \code{\link{query_relationship}}
#' 
#' @examples
#' tped <- tidyped(small_ped)
#' 
#' # Compact matrix
#' A_compact <- pedmat(tped, method = "A", compact = TRUE)
#' dim(A_compact)  # Reduced dimensions
#' 
#' # Expand to full size
#' A_full <- expand_pedmat(A_compact)
#' dim(A_full)  # Original dimensions restored
#' 
#' # Non-compact matrices are returned unchanged
#' A <- pedmat(tped, method = "A", compact = FALSE)
#' A2 <- expand_pedmat(A)
#' identical(dim(A), dim(A2))  # TRUE
#' 
#' @export
expand_pedmat <- function(x) {
  # NOTE: sibling off-diagonal correction uses calc_sib_offdiag() below.
  # Check if x is a pedmat (either S3 or S4 with marker)
  is_pedmat <- inherits(x, "pedmat") || !is.null(attr(x, "pedmat_S4"))
  if (!is_pedmat) {
    stop("x must be a pedmat object")
  }
  
  ci <- attr(x, "call_info")
  if (!isTRUE(ci$compact)) {
    # Already full size, just return the matrix
    result <- x
    if (inherits(result, "pedmat")) {
      class(result) <- setdiff(class(result), "pedmat")
    }
    return(result)
  }
  
  map <- attr(x, "compact_map")
  mat <- x
  if (inherits(mat, "pedmat")) {
    class(mat) <- setdiff(class(mat), "pedmat")
  }
  
  # Get original order based on Ind column
  # compact_map is sorted by OldIndNum (original pedigree order)
  # We need to get the original pedigree Ind order
  ped_full <- attr(x, "call_info")$ped_original
  if (is.null(ped_full)) {
    # Fallback: use map$Ind order as-is
    orig_ind_order <- map$Ind
  } else {
    orig_ind_order <- ped_full$Ind
  }
  
  # Filter map to match original order
  map_ordered <- map[match(orig_ind_order, map$Ind)]
  
  # Matrix-like or vector-like?
  if (!is.null(dim(mat))) {
    # Matrix expansion
    result_full <- mat[map_ordered$RepIndNum, map_ordered$RepIndNum]
    rownames(result_full) <- map_ordered$Ind
    colnames(result_full) <- map_ordered$Ind
    
    # Fix sibling off-diagonals for A, D, AA
    # In compact mode, full-siblings share a representative. Their relationship to
    # others is identical, but their relationship to each other (off-diagonal) 
    # is different from their relationship to self (diagonal).
    primary_method <- ci$method[1]
    if (primary_method %in% c("A", "D", "AA")) {
      counts <- table(map_ordered$RepIndNum)
      reps_with_sibs <- as.integer(names(counts[counts > 1]))
      
      if (length(reps_with_sibs) > 0) {
        # Get A matrix (required for sibling calculation)
        A_mat <- if (primary_method == "A") mat else attr(x, "A_intermediate")
        
        if (!is.null(A_mat)) {
          for (rep_idx in reps_with_sibs) {
            first_idx <- which(map_ordered$RepIndNum == rep_idx)[1]
            s_idx <- map_ordered$SireNum[first_idx]
            d_idx <- map_ordered$DamNum[first_idx]
            
            sib_val <- calc_sib_offdiag(A_mat, s_idx, d_idx, primary_method)
            
            member_indices <- which(map_ordered$RepIndNum == rep_idx)
            if (length(member_indices) > 1) {
              diag_val <- mat[rep_idx, rep_idx]
              result_full[member_indices, member_indices] <- sib_val
              diag_indices <- cbind(member_indices, member_indices)
              result_full[diag_indices] <- diag_val
            }
          }
        }
      }
    }
    return(result_full)
  } else {
    # Vector expansion
    result_full <- mat[map_ordered$RepIndNum]
    names(result_full) <- map_ordered$Ind
    return(result_full)
  }
}

# --- S3 Methods for pedmat class ---

#' @export
print.pedmat <- function(x, ...) {
  ci <- attr(x, "call_info")
  cat("Pedigree Matrix (", paste(ci$method, collapse=", "), ")\n", sep="")
  if (isTRUE(ci$compact)) {
    cat("Compacted: ", ci$n_original, " -> ", ci$n_compact, " individuals\n", sep="")
  } else {
    cat("Size: ", ci$n_original, " individuals\n", sep="")
  }
  
  cat("\nUse summary() for details, attr(x, 'ped') for pedigree data.\n\n")
  
  # Strip pedmat class and print using underlying method
  x_show <- x
  if (inherits(x_show, "pedmat")) {
    class(x_show) <- setdiff(class(x_show), "pedmat")
  }
  print(x_show, ...)
}

#' Summary Statistics for Pedigree Matrices
#' 
#' @description
#' Computes and displays summary statistics for a pedmat object.
#' 
#' @param x A pedmat object from \code{\link{pedmat}}.
#' 
#' @return An object of class \code{"summary.pedmat"} with statistics including
#'   method, dimensions, compression ratio (if compact), mean relationship, 
#'   and matrix density.
#' 
#' @details
#' Since pedmat objects are often S4 sparse matrices with custom attributes,
#' use this function instead of the generic \code{summary()} to ensure proper
#' display of pedigree matrix statistics.
#' 
#' @examples
#' tped <- tidyped(small_ped)
#' A <- pedmat(tped, method = "A")
#' summary_pedmat(A)
#' 
#' @seealso \code{\link{pedmat}}
#' @export
summary_pedmat <- function(x) {
  # Check if it's a pedmat
  is_pedmat <- inherits(x, "pedmat") || !is.null(attr(x, "pedmat_S4"))
  if (!is_pedmat) {
    stop("x must be a pedmat object")
  }
  
  ci <- attr(x, "call_info")
  
  stats <- list(
    method = ci$method,
    compact = isTRUE(ci$compact),
    n_original = ci$n_original,
    n_calculated = ci$n_compact
  )
  
  # Calculate matrix-level stats
  obj_clean <- x
  class(obj_clean) <- setdiff(class(obj_clean), "pedmat")
  
  if (!is.list(obj_clean)) {
    # For Matrix objects, use Matrix::mean or extract values
    stats$mean_val <- tryCatch({
      n <- as.numeric(nrow(obj_clean))
      denom <- n * n - n
      if (denom <= 0) return(NA_real_)
      if (inherits(obj_clean, "Matrix")) {
        total_sum <- sum(obj_clean)
        diag_sum <- sum(Matrix::diag(obj_clean))
        (total_sum - diag_sum) / denom
      } else {
        total_sum <- sum(obj_clean)
        diag_sum <- sum(diag(obj_clean))
        (total_sum - diag_sum) / denom
      }
    }, error = function(e) NA_real_)
    
    stats$sparsity <- if (inherits(obj_clean, "Matrix")) {
      # Works for both sparseMatrix (e.g. dgCMatrix) and dense Matrix (e.g. dgeMatrix).
      # speed branch returns dgeMatrix for A/D/AA to avoid the N² zero-scan
      # overhead of Matrix::Matrix(..., sparse=TRUE); nnzero() handles both cases.
      Matrix::nnzero(obj_clean) / (as.numeric(nrow(obj_clean)) * ncol(obj_clean))
    } else if (is.matrix(obj_clean)) {
      sum(obj_clean != 0L) / length(obj_clean)
    } else { NA_real_ }
  }
  
  if (isTRUE(ci$compact)) {
    stats$compact_stats <- attr(x, "compact_stats")
    fs <- attr(x, "family_summary")
    if (!is.null(fs)) {
      stats$top_families <- head(fs[order(-fs$FamilySize), ], 5)
    }
  }
  
  class(stats) <- "summary.pedmat"
  return(stats)
}

#' @export
summary.pedmat <- function(object, ...) {
  # Delegate to summary_pedmat
  summary_pedmat(object)
}

#' @export
print.summary.pedmat <- function(x, ...) {
  cat("Summary of Pedigree Matrix (", paste(x$method, collapse=", "), ")\n", sep="")
  cat("========================================\n")
  cat("Input Size:     ", x$n_original, " individuals\n")
  cat("Calculated Size:", x$n_calculated, " individuals\n")
  
  if (x$compact) {
    cat("\nCompaction Results:\n")
    if (!is.null(x$compact_stats)) {
      ratio <- x$compact_stats$compression_ratio
      cat("- Compression: ", round(ratio * 100, 1), "%\n", sep="")
      cat("- Families:    ", x$compact_stats$n_families_compacted, " families merged\n")
    }
    if (!is.null(x$top_families)) {
      cat("\nTop 5 families by size:\n")
      print(x$top_families[, c("FamilyLabel", "FamilySize")])
    }
  }
  
  if (!is.null(x$mean_val)) {
      cat("\nMatrix Properties:\n")
      cat("- Mean off-diagonal relationship: ", if(is.na(x$mean_val)) "NA" else round(x$mean_val, 6), "\n")
      cat("- Density (non-zero): ", if(is.na(x$sparsity)) "NA" else round(x$sparsity * 100, 2), "%\n", sep="")
  }
  cat("========================================\n")
}

# --- Generic functions to enable S3 dispatch ---



#' @export
#' @method [ pedmat
`[.pedmat` <- function(x, ...) {
  # Strip pedmat class and delegate to underlying method
  class(x) <- setdiff(class(x), "pedmat")
  x[...]
}

#' @export
#' @method dim pedmat
dim.pedmat <- function(x) {
  class(x) <- setdiff(class(x), "pedmat")
  dim(x)
}

#' @export
#' @method length pedmat
length.pedmat <- function(x) {
  class(x) <- setdiff(class(x), "pedmat")
  length(x)
}

#' @export
#' @method t pedmat
t.pedmat <- function(x) {
  class(x) <- setdiff(class(x), "pedmat")
  t(x)
}

#' @export
as.matrix.pedmat <- function(x, ...) {
  class(x) <- setdiff(class(x), "pedmat")
  if (inherits(x, "Matrix")) return(as.matrix(x))
  return(x)
}

#' @export
as.vector.pedmat <- function(x, mode = "any") {
  class(x) <- setdiff(class(x), "pedmat")
  as.vector(x, mode = mode)
}
