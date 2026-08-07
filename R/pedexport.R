# ---- pedexport: Export tidyped pedigrees to external breeding software formats ----

#' Export a pedigree to a breeding software format
#'
#' Converts a \code{tidyped} pedigree into the tabular format expected by
#' common animal/plant breeding software packages and either returns the
#' result as a \code{data.table} or writes it to a plain-text file.
#'
#' @param ped A complete \code{tidyped} object (see \code{\link{tidyped}}).
#'   The pedigree must be structurally complete: every non-missing \code{Sire}
#'   and \code{Dam} identifier must be present in \code{Ind}.
#' @param software Character scalar specifying the target format.  One of:
#'   \describe{
#'     \item{\code{"blupf90"}}{Three integer columns (\code{IndNum},
#'       \code{SireNum}, \code{DamNum}), no header, missing parents encoded
#'       as \code{0}.  Rows are sorted by \code{IndNum} so parent rows always
#'       precede offspring rows.  Compatible with BLUPF90, REMLF90,
#'       GIBBSF90, and related programs.}
#'     \item{\code{"asreml"}}{Three character columns (\code{animal},
#'       \code{sire}, \code{dam}), with a one-line header, missing parents
#'       encoded as \code{"0"}.  Rows are sorted by generation (\code{Gen})
#'       so founders appear first.  Compatible with ASReml-R and ASReml-SA.}
#'     \item{\code{"wombat"}}{Identical integer layout to \code{"blupf90"}.
#'       Compatible with Wombat.}
#'     \item{\code{"mtdfreml"}}{Identical integer layout to \code{"blupf90"}.
#'       Compatible with MTDFREML and MTGSAM.}
#'     \item{\code{"numeric"}}{Integer layout with an optional header
#'       (\code{IndNum SireNum DamNum}).  A portable, software-agnostic
#'       choice when the target program is not listed above.}
#'   }
#' @param file Character scalar.  Path to the output file.  When \code{NULL}
#'   (default) the formatted \code{data.table} is returned invisibly and no
#'   file is written.
#' @param sep Character scalar.  Field separator used when writing to
#'   \code{file}.  Defaults to a single space (\code{" "}).
#' @param header Logical scalar or \code{NULL}.  Whether to include a column
#'   header line.  \code{NULL} (default) uses the software-specific default:
#'   \code{TRUE} for \code{"asreml"} and \code{"numeric"}, \code{FALSE} for
#'   all others.
#' @param missing Character or integer scalar.  Symbol for missing parents.
#'   \code{NULL} (default) uses the software-specific default: \code{0L} for
#'   numeric formats, \code{"0"} for \code{"asreml"}.
#'
#' @return A \code{data.table} in the target format, returned invisibly.
#'   When \code{file} is not \code{NULL}, the table is also written to that
#'   path and a message reports the number of individuals written.
#'
#' @details
#' \strong{Column semantics by format:}
#'
#' \tabular{llll}{
#'   \strong{software} \tab \strong{Col 1} \tab \strong{Col 2} \tab
#'     \strong{Col 3} \cr
#'   blupf90 / wombat / mtdfreml / numeric \tab \code{IndNum} (integer) \tab
#'     \code{SireNum} (integer) \tab \code{DamNum} (integer) \cr
#'   asreml \tab \code{animal} (character) \tab \code{sire} (character) \tab
#'     \code{dam} (character) \cr
#' }
#'
#' \strong{Row order:}
#'
#' All numeric formats sort by \code{IndNum} ascending, which guarantees that
#' parent rows appear before offspring rows (parents always receive a smaller
#' integer index after \code{tidyped()} topological sorting).  The
#' \code{"asreml"} format sorts by \code{Gen} ascending for the same reason.
#'
#' \strong{Relationship to \code{tidyped()}:}
#'
#' \code{pedexport()} is a post-processing step that re-encodes an already
#' validated pedigree.  Pass the output of \code{\link{tidyped}} directly to
#' avoid redundant validation.
#'
#' @seealso \code{\link{tidyped}} to prepare the pedigree,
#'   \code{\link{pedmat}} to compute relationship matrices.
#'
#' @examples
#' library(visPedigree)
#'
#' tp <- tidyped(small_ped)
#'
#' # Return a data.table (no file written)
#' out_blupf90 <- pedexport(tp, software = "blupf90")
#' head(out_blupf90)
#'
#' # ASReml format with character IDs and header
#' out_asreml <- pedexport(tp, software = "asreml")
#' head(out_asreml)
#'
#' \donttest{
#' # Write to a file
#' tmp <- tempfile(fileext = ".txt")
#' pedexport(tp, software = "blupf90", file = tmp)
#' readLines(tmp, n = 5)
#' }
#'
#' @import data.table
#' @export
pedexport <- function(ped,
                      software = c("blupf90", "asreml", "wombat",
                                   "mtdfreml", "numeric"),
                      file    = NULL,
                      sep     = " ",
                      header  = NULL,
                      missing = NULL) {

  # ---- 1. Input validation ----
  ped      <- ensure_complete_tidyped(ped, "pedexport()")
  software <- match.arg(software)

  if (!is.null(file)) {
    if (!is.character(file) || length(file) != 1L || nchar(file) == 0L)
      stop("'file' must be a non-empty character string or NULL.", call. = FALSE)
  }

  if (!is.character(sep) || length(sep) != 1L)
    stop("'sep' must be a single character string.", call. = FALSE)

  if (!is.null(header) && (!is.logical(header) || length(header) != 1L || is.na(header)))
    stop("'header' must be TRUE, FALSE, or NULL.", call. = FALSE)

  # ---- 2. Software-specific defaults ----
  use_char   <- software == "asreml"
  def_header <- use_char || software == "numeric"
  def_miss   <- if (use_char) "0" else 0L

  if (is.null(header))  header  <- def_header
  if (is.null(missing)) missing <- def_miss

  # ---- 3. Build output table ----
  out <- if (use_char) {
    .pedexport_char(ped, missing = as.character(missing))
  } else {
    .pedexport_num(ped, missing = as.integer(missing))
  }

  # ---- 4. Write file if requested ----
  if (!is.null(file)) {
    data.table::fwrite(out, file = file, sep = sep, col.names = header,
                       quote = FALSE)
    message(sprintf("Written %d individuals to: %s", nrow(out), file))
  }

  invisible(out[])
}


# ---- Internal helpers --------------------------------------------------------

#' Build a numeric (integer) export table
#' @param ped A complete tidyped object.
#' @param missing Integer scalar for missing parents.
#' @return A data.table with columns IndNum, SireNum, DamNum.
#' @keywords internal
.pedexport_num <- function(ped, missing = 0L) {
  out <- data.table::data.table(
    IndNum  = ped$IndNum,
    SireNum = ped$SireNum,
    DamNum  = ped$DamNum
  )
  # Replace 0 (internal missing sentinel) with the requested symbol
  if (missing != 0L) {
    out[SireNum == 0L, SireNum := missing]
    out[DamNum  == 0L, DamNum  := missing]
  }
  data.table::setorderv(out, "IndNum")
  out[]
}

#' Build a character export table (ASReml format)
#' @param ped A complete tidyped object.
#' @param missing Character scalar for missing parents.
#' @return A data.table with columns animal, sire, dam.
#' @keywords internal
.pedexport_char <- function(ped, missing = "0") {
  sire_char <- ped$Sire
  dam_char  <- ped$Dam
  sire_char[is.na(sire_char)] <- missing
  dam_char[is.na(dam_char)]   <- missing

  out <- data.table::data.table(
    animal = ped$Ind,
    sire   = sire_char,
    dam    = dam_char,
    Gen    = ped$Gen
  )
  data.table::setorderv(out, "Gen")
  out[, Gen := NULL]
  out[]
}
