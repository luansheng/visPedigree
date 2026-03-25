#' Calculate inbreeding coefficients
#'
#' \code{inbreed} function calculates the inbreeding coefficients for all individuals in a tidied pedigree.
#'
#' This function takes a pedigree tidied by the \code{\link{tidyped}} function and
#' calculates the inbreeding coefficients using an optimized C++ implementation of
#' the Sargolzaei & Iwaisaki (2005) LAP (Longest Ancestral Path) bucket algorithm.
#' This method is the fastest known direct algorithm for computing all inbreeding
#' coefficients: it replaces the O(\eqn{N^2}) linear scan of Meuwissen & Luo (1992)
#' with O(1) bucket pops and selective ancestor clearing, giving
#' \eqn{O(\sum m_i)} total work where \eqn{m_i} is the number of distinct
#' ancestors of individual \eqn{i}. At \eqn{N = 1{,}000{,}000}, the kernel
#' completes in approximately 0.12 s — over 10\eqn{\times} faster than the previous
#' Meuwissen & Luo (1992) implementation and on par with the \pkg{pedigreemm}
#' reference C implementation of the same algorithm.
#' It is the core engine used by both \code{tidyped(..., inbreed = TRUE)} and
#' \code{pedmat(..., method = "f")}, ensuring consistent results across the package.
#'
#' @param ped A \code{tidyped} object.
#' @param ... Additional arguments (currently ignored).
#' @return A \code{tidyped} object with an additional column \strong{f}.
#' @examples
#' library(visPedigree)
#' data(simple_ped)
#' ped <- tidyped(simple_ped)
#' ped_f <- inbreed(ped)
#' ped_f[f > 0, .(Ind, Sire, Dam, f)]
#' @export
#' @import data.table
inbreed <- function(ped, ...) {
  ped <- ensure_complete_tidyped(ped, "inbreed()")

  meta <- attr(ped, "ped_meta")

  # tidyped guarantees IndNum == seq_len(N) and rows are sorted by IndNum,
  # so no re-sorting is needed before calling the C++ engine.
  # We still copy to avoid mutating the caller's object in place.
  ped_work <- copy(ped)

  res_f <- cpp_calculate_inbreeding(ped_work$SireNum, ped_work$DamNum)

  ped_work[, f := res_f$f]

  # Propagate ped_meta from original object
  result <- new_tidyped(ped_work)
  if (!is.null(meta)) data.table::setattr(result, "ped_meta", meta)

  return(result)
}
