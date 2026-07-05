#' Matrix-Free Pedigree Relationship Products
#'
#' Computes products of the form \eqn{A x}, \eqn{A X}, \eqn{A^{-1} x}, or
#' \eqn{A^{-1} X} directly from an ordered pedigree, without constructing the
#' full additive relationship matrix \eqn{A}. This is useful for large
#' pedigrees where explicitly storing the dense \eqn{n \times n} matrix is
#' infeasible.
#'
#' @param ped A complete pedigree accepted by \code{\link{tidyped}}, preferably
#'   an existing \code{tidyped} object.
#' @param x A numeric or logical vector, or a numeric or logical matrix.
#'   An unnamed vector must have one value per pedigree individual in
#'   \code{IndNum} order. A named vector may contain a subset of individuals;
#'   omitted individuals receive value zero. For a matrix, the same rules apply
#'   to its rows and \code{rownames}.
#' @param method Character scalar. \code{"A"} (default) computes products with
#'   the additive relationship matrix. \code{"Ainv"} computes products with its
#'   inverse.
#'
#' @details
#' For an ordered pedigree, the additive relationship matrix can be written as
#' \deqn{A = T D T^\prime,}
#' where \eqn{T} is the transmission matrix and \eqn{D} contains Mendelian
#' sampling variances. \code{pedprod()} applies these factors through backward
#' and forward pedigree traversals. For an \eqn{n \times p} right-hand side,
#' computation requires \eqn{O(np)} time and working memory rather than
#' materializing the \eqn{O(n^2)} relationship matrix.
#'
#' Named inputs are aligned by individual ID. Every supplied name must occur in
#' \code{ped$Ind}, and duplicate, missing, or empty names are rejected. Missing
#' values and non-finite input values are also rejected to avoid propagating an
#' invalid value through all descendants or relatives.
#'
#' @return For vector input, a named numeric vector in pedigree order. For
#'   matrix input, a numeric matrix whose rows are in pedigree order; input
#'   column names are preserved.
#'
#' @references
#' Colleau, J. J. (2002). An indirect approach to the extensive calculation of
#' relationship coefficients. \emph{Genetics Selection Evolution}, 34,
#' 409--421. \doi{10.1051/gse:2002015}
#'
#' Colleau, J. J., Palhiere, I., Rodriguez-Ramilo, S. T., & Legarra, A.
#' (2017). A fast indirect method to compute functions of genomic relationships
#' concerning genotyped and ungenotyped individuals, for diversity management.
#' \emph{Genetics Selection Evolution}, 49, 87.
#' \doi{10.1186/s12711-017-0363-9}
#'
#' @seealso \code{\link{pedmat}} for explicitly constructing relationship
#'   matrices and \code{\link{pedrel}} for group relationship summaries.
#'
#' @examples
#' tped <- tidyped(small_ped)
#'
#' # Equal contributions from two candidates; unspecified individuals are zero.
#' contributions <- c(Z1 = 0.5, Z2 = 0.5)
#' relationship_to_group <- pedprod(tped, contributions)
#'
#' # Average additive relationship among the weighted candidates: c' A c.
#' sum(contributions * relationship_to_group[names(contributions)])
#'
#' # Multiple contribution schemes can be evaluated in one call.
#' schemes <- cbind(
#'   Equal = c(A = 0.5, B = 0.5),
#'   SireA = c(A = 1.0, B = 0.0)
#' )
#' pedprod(tped, schemes)
#'
#' @export
pedprod <- function(ped, x, method = c("A", "Ainv")) {
  method <- match.arg(method)

  if (inherits(ped, "splitped")) {
    stop(
      "pedprod() does not support 'splitped' objects directly. ",
      "Apply pedprod() to one complete pedigree group at a time.",
      call. = FALSE
    )
  }

  ped <- data.table::copy(ped)
  if (is_tidyped(ped) ||
      all(c("Ind", "Sire", "Dam", "Sex", "Gen",
            "IndNum", "SireNum", "DamNum") %in% names(ped))) {
    ped <- ensure_complete_tidyped(ped, "pedprod()")
  } else {
    ped <- tidyped(ped, addnum = TRUE)
  }

  if (!all(c("IndNum", "SireNum", "DamNum") %in% names(ped))) {
    ped[, IndNum := .I]
    ped[, SireNum := match(Sire, Ind, nomatch = 0L)]
    ped[, DamNum := match(Dam, Ind, nomatch = 0L)]
  }

  data.table::setorder(ped, IndNum)
  rhs <- prepare_pedprod_rhs(x, ped$Ind)
  f_res <- cpp_calculate_inbreeding(ped$SireNum, ped$DamNum)

  result <- switch(
    method,
    A = cpp_multiply_A(ped$SireNum, ped$DamNum, f_res$dii, rhs$values),
    Ainv = cpp_multiply_Ainv(ped$SireNum, ped$DamNum, f_res$dii, rhs$values)
  )

  rownames(result) <- ped$Ind
  colnames(result) <- rhs$colnames

  if (rhs$vector) {
    result <- as.numeric(result[, 1L])
    names(result) <- ped$Ind
  }

  result
}

# Prepare a vector or matrix right-hand side in pedigree row order.
prepare_pedprod_rhs <- function(x, ids) {
  is_vector <- is.atomic(x) && is.null(dim(x))
  is_matrix <- is.matrix(x)

  if (!is_vector && !is_matrix) {
    stop("'x' must be a numeric or logical vector or matrix.", call. = FALSE)
  }
  if (!is.numeric(x) && !is.logical(x)) {
    stop("'x' must be numeric or logical.", call. = FALSE)
  }
  if (length(x) == 0L || (is_matrix && ncol(x) == 0L)) {
    stop("'x' must contain at least one value.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop("'x' must contain only finite, non-missing values.", call. = FALSE)
  }

  n <- length(ids)

  if (is_vector) {
    x_names <- names(x)
    if (is.null(x_names)) {
      if (length(x) != n) {
        stop(
          sprintf("Unnamed 'x' must have length %d; got %d.", n, length(x)),
          call. = FALSE
        )
      }
      values <- matrix(as.numeric(x), ncol = 1L)
    } else {
      validate_pedprod_ids(x_names, ids, "names(x)")
      values <- matrix(0, nrow = n, ncol = 1L)
      values[match(x_names, ids), 1L] <- as.numeric(x)
    }

    return(list(values = values, vector = TRUE, colnames = NULL))
  }

  x_rows <- rownames(x)
  if (is.null(x_rows)) {
    if (nrow(x) != n) {
      stop(
        sprintf("Unlabelled matrix 'x' must have %d rows; got %d.", n, nrow(x)),
        call. = FALSE
      )
    }
    values <- matrix(
      as.numeric(x),
      nrow = nrow(x),
      ncol = ncol(x),
      dimnames = list(NULL, colnames(x))
    )
  } else {
    validate_pedprod_ids(x_rows, ids, "rownames(x)")
    values <- matrix(
      0,
      nrow = n,
      ncol = ncol(x),
      dimnames = list(NULL, colnames(x))
    )
    values[match(x_rows, ids), ] <- x
  }

  list(values = values, vector = FALSE, colnames = colnames(x))
}

# Validate user-supplied individual identifiers before zero-filling/alignment.
validate_pedprod_ids <- function(input_ids, pedigree_ids, label) {
  if (anyNA(input_ids) || any(input_ids == "")) {
    stop(sprintf("%s must not contain missing or empty IDs.", label), call. = FALSE)
  }
  if (anyDuplicated(input_ids)) {
    stop(sprintf("%s must contain unique individual IDs.", label), call. = FALSE)
  }

  unknown <- setdiff(input_ids, pedigree_ids)
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "%s contains unknown individual IDs: %s.",
        label,
        paste(utils::head(unknown, 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
