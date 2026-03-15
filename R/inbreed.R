#' Calculate inbreeding coefficients
#'
#' \code{inbreed} function calculates the inbreeding coefficients for all individuals in a tidied pedigree.
#'
#' This function takes a pedigree tidied by the \code{\link{tidyped}} function and calculates the inbreeding coefficients using optimized C++ code based on the Meuwissen & Luo (1992) algorithm. It is the core engine used by both \code{tidyped(..., inbreed = TRUE)} and \code{pedmat(..., method = "f")}, ensuring consistent results across the package. It is significantly faster than standard R implementations for large pedigrees.
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
  ped_work <- copy(ped)
  ped_work[, .orig_order_visped := seq_len(.N)]
  setorder(ped_work, IndNum)

  # Call Rcpp function
  res_f <- cpp_calculate_inbreeding(ped_work$SireNum, ped_work$DamNum)
  
  ped_work[, f := res_f$f]
  setorder(ped_work, .orig_order_visped)
  ped_work[, .orig_order_visped := NULL]
  
  # Propagate ped_meta from original object
  result <- new_tidyped(ped_work)
  if (!is.null(meta)) data.table::setattr(result, "ped_meta", meta)
  
  return(result)
}
