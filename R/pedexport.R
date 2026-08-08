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
#'       \code{sire}, \code{dam}), with a one-line header by default, missing
#'       parents encoded as \code{"0"}.  Rows are sorted by generation
#'       (\code{Gen}) so founders appear first.  Compatible with ASReml-R
#'       and ASReml-SA (character IDs require the \code{!ALPHA} qualifier;
#'       the header line requires the \code{!SKIP 1} qualifier).}
#'     \item{\code{"wombat"}}{Three character columns (\code{animal},
#'       \code{sire}, \code{dam}), no header, missing parents encoded as
#'       \code{"0"}.  Wombat accepts alphanumeric IDs and recodes them
#'       internally.  Rows are sorted so parents precede offspring.}
#'     \item{\code{"mtdfreml"}}{Identical integer layout to \code{"blupf90"}.
#'       Compatible with MTDFREML and MTGSAM.}
#'     \item{\code{"dmu"}}{Identical integer layout to \code{"blupf90"}.
#'       Compatible with DMU.}
#'     \item{\code{"numeric"}}{Integer layout with an optional header
#'       (\code{IndNum SireNum DamNum}).  A portable, software-agnostic
#'       choice when the target program is not listed above.}
#'     \item{\code{"sommer"}}{Three character columns (\code{ID}, \code{Sire},
#'       \code{Dam}), missing parents coded as \code{NA}.  Rows are sorted so
#'       parents precede offspring.  Returned as a data.table ready to pair
#'       with \code{\link{pedmat}} for the relationship matrix used by
#'       \code{sommer::mmer(..., random = ~vsr(animal, Gu = A))}.}
#'   }
#' @param file Character scalar.  Path to the output file.  When \code{NULL}
#'   (default) the formatted \code{data.table} is returned invisibly and no
#'   file is written.  File output is not supported for \code{"sommer"},
#'   which is intended for direct use as an R \code{data.table}.
#' @param sep Character scalar.  Field separator used when writing to
#'   \code{file}.  Defaults to a single space (\code{" "}), which every
#'   file-based format accepts.  BLUPF90 requires a space; WOMBAT, MTDFREML,
#'   and DMU accept spaces or TABs; ASReml accepts spaces, TABs, or commas;
#'   and the generic \code{"numeric"} format accepts any single-byte,
#'   non-empty separator supported by \code{\link[data.table]{fwrite}}.
#'   ASReml comma-delimited output requires a \code{.csv} file extension or
#'   the \code{!CSV} qualifier.
#' @param header Logical scalar or \code{NULL}.  Whether to include a column
#'   header line.  \code{NULL} (default) uses the software-specific default:
#'   \code{TRUE} for \code{"asreml"} and \code{"numeric"},
#'   \code{FALSE} for \code{"blupf90"}, \code{"wombat"}, \code{"mtdfreml"}
#'   and \code{"dmu"}.  Ignored for \code{"sommer"}.  Note for ASReml: a
#'   header line is read as data unless the \code{!SKIP 1} qualifier is used
#'   in the ASReml command file.
#' @param missing Character or integer scalar.  Symbol for missing parents.
#'   \code{NULL} (default) uses the software-specific default: \code{0L} for
#'   numeric formats, \code{"0"} for \code{"asreml"} and \code{"wombat"}.
#'   Numeric formats (\code{"blupf90"}, \code{"mtdfreml"}, \code{"dmu"},
#'   \code{"numeric"}) require a single integer value; \code{"asreml"} and
#'   \code{"wombat"} accept a character value (numeric values are converted
#'   to character).  Ignored for \code{"sommer"}, which always codes missing
#'   parents as \code{NA}.
#'
#' @return A \code{data.table} in the target format, returned invisibly.
#'   Numeric formats carry an \code{xref} attribute mapping each numeric ID
#'   back to its original character ID (columns \code{IndNum} and
#'   \code{Ind}), and when \code{file} is given the mapping is also written
#'   to \code{paste0(file, ".xref")} without a header, mirroring RENUMF90's
#'   \code{_XrefID} file.  For file-based formats, when \code{file} is not
#'   \code{NULL}, the table is also written to that path and a message reports
#'   the number of individuals written.
#'
#' @details
#' \strong{Column semantics by format:}
#'
#' \tabular{llll}{
#'   \strong{software} \tab \strong{Col 1} \tab \strong{Col 2} \tab
#'     \strong{Col 3} \cr
#'   blupf90 / mtdfreml / dmu / numeric \tab \code{IndNum} (integer) \tab
#'     \code{SireNum} (integer) \tab \code{DamNum} (integer) \cr
#'   asreml / wombat \tab \code{animal} (character) \tab \code{sire} (character) \tab
#'     \code{dam} (character) \cr
#'   sommer \tab \code{ID} (character) \tab \code{Sire} (character) \tab
#'     \code{Dam} (character) \cr
#' }
#'
#' \strong{Software format requirements:}
#'
#' All file-based formats use unquoted, free-format plain text.  The table
#' below summarises what each program expects, and the defaults
#' \code{pedexport()} uses to meet those expectations.
#'
#' \tabular{llllll}{
#'   \strong{software} \tab \strong{header} \tab \strong{separator} \tab
#'     \strong{missing} \tab \strong{IDs} \tab \strong{notes} \cr
#'   blupf90 \tab no \tab spaces only \tab \code{0} \tab integer \tab
#'     TAB separators rejected \cr
#'   wombat \tab no \tab space / TAB \tab \code{"0"} \tab character \tab
#'     alphanumeric accepted, recoded internally \cr
#'   mtdfreml \tab no \tab space / TAB \tab \code{0} \tab integer \tab \cr
#'   dmu \tab no \tab space / TAB \tab \code{0} \tab integer \tab \cr
#'   asreml \tab yes \tab space / TAB / comma \tab \code{"0"} \tab
#'     char or integer \tab character IDs need \code{!ALPHA}; header needs
#'     \code{!SKIP 1}; comma needs \code{.csv} or \code{!CSV} \cr
#'   sommer \tab n/a \tab n/a \tab \code{NA} \tab character \tab returned
#'     as a data.table, not a file \cr
#'   numeric \tab optional \tab any \tab \code{0} \tab integer \tab generic
#'     software-agnostic layout \cr
#' }
#'
#' \strong{Row order:}
#'
#' All numeric formats sort by \code{IndNum} ascending, which guarantees that
#' parent rows appear before offspring rows (parents always receive a smaller
#' integer index after \code{tidyped()} topological sorting).  The
#' \code{"asreml"}, \code{"wombat"} and \code{"sommer"} formats sort by
#' \code{Gen} ascending for the same reason.
#'
#' \strong{Optional \code{tidyped()} columns:}
#'
#' \code{tidyped(..., addnum = FALSE)} omits \code{IndNum}/\code{SireNum}/
#' \code{DamNum} and \code{tidyped(..., addgen = FALSE)} omits \code{Gen}.
#' When these columns are absent, \code{pedexport()} reconstructs them from
#' \code{Ind}/\code{Sire}/\code{Dam}.  Because \code{tidyped()} always
#' topologically sorts rows, parents still precede offspring in the output.
#'
#' \strong{Identifiers in written files:}
#'
#' Output is deliberately unquoted because the supported breeding programs
#' read pedigree files as free-format text.  Therefore, identifiers and
#' character missing-parent symbols must not contain whitespace.  They also
#' must not contain the selected non-whitespace separator.  Numeric formats
#' apply the same restriction to original identifiers written to the
#' companion \code{.xref} file.
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
#' # Numeric formats carry the ID mapping back to character IDs
#' out_dmu <- pedexport(tp, software = "dmu")
#' head(attr(out_dmu, "xref"))
#'
#' # sommer: character pedigree + visPedigree's own relationship matrix.
#' # pedmat() rownames match pedexport(..., software = "sommer")$ID, so the
#' # matrix is ready for sommer::mmer(y ~ 1, random = ~vsr(animal, Gu = A)).
#' out_sommer <- pedexport(tp, software = "sommer")
#' head(out_sommer)
#' A <- pedmat(tp, method = "A")
#'
#' \donttest{
#' # Write to a file (numeric formats also write <file>.xref)
#' tmp <- tempfile(fileext = ".txt")
#' pedexport(tp, software = "blupf90", file = tmp)
#' readLines(tmp, n = 5)
#' readLines(paste0(tmp, ".xref"), n = 5)
#' }
#'
#' @import data.table
#' @export
pedexport <- function(ped,
                      software = c("blupf90", "asreml", "wombat",
                                   "mtdfreml", "dmu", "numeric", "sommer"),
                      file    = NULL,
                      sep     = " ",
                      header  = NULL,
                      missing = NULL) {

  # ---- 1. Input validation ----
  ped      <- ensure_complete_tidyped(ped, "pedexport()")
  software <- match.arg(software)

  if (!is.null(file)) {
    if (!is.character(file) || length(file) != 1L || is.na(file) ||
        !nzchar(file)) {
      stop("'file' must be a non-empty character string or NULL.", call. = FALSE)
    }
    if (software == "sommer") {
      stop("File output is not supported for software = \"sommer\"; ",
           "use the returned data.table directly.", call. = FALSE)
    }
  }

  if (!is.character(sep) || length(sep) != 1L || is.na(sep) ||
      nchar(sep, type = "bytes") != 1L || sep %in% c("\r", "\n")) {
    stop("'sep' must be one non-empty, non-newline, single-byte character.",
         call. = FALSE)
  }
  .validate_export_sep(software, sep)

  if (!is.null(header) && (!is.logical(header) || length(header) != 1L || is.na(header)))
    stop("'header' must be TRUE, FALSE, or NULL.", call. = FALSE)

  # ---- 2. Software-specific defaults ----
  use_char   <- software %in% c("asreml", "wombat", "sommer")
  def_header <- software %in% c("asreml", "numeric")
  def_miss   <- if (use_char) "0" else 0L

  if (is.null(header))  header  <- def_header
  if (is.null(missing)) missing <- def_miss

  # ---- 3. Validate 'missing' per format ----
  if (software == "sommer") {
    # R convention: missing parents are always NA; the symbol is not used.
    missing <- NA_character_
  } else if (use_char) {
    if ((!is.character(missing) && !is.numeric(missing)) ||
        length(missing) != 1L || is.na(missing)) {
      stop("For the 'asreml' and 'wombat' formats, 'missing' must be a ",
           "single character (or numeric) value.", call. = FALSE)
    }
    missing <- as.character(missing)
    if (!nzchar(missing)) {
      stop("For the 'asreml' and 'wombat' formats, 'missing' must not be empty.",
           call. = FALSE)
    }
  } else {
    missing_int <- suppressWarnings(as.integer(missing))
    if (!is.numeric(missing) || length(missing) != 1L || is.na(missing) ||
        is.na(missing_int) || missing_int != missing) {
      stop("For numeric formats (blupf90, mtdfreml, dmu, numeric), ",
           "'missing' must be a single integer value.", call. = FALSE)
    }
    missing <- missing_int
  }

  # ---- 4. Build output table ----
  # tidyped(..., addnum = FALSE) / addgen = FALSE omit the ordering columns;
  # reconstruct them from Ind/Sire/Dam (see .ensure_export_cols).
  ped <- .ensure_export_cols(ped)
  out <- if (use_char) {
    .pedexport_char(ped, missing = missing,
                    cols = if (software == "sommer") c("ID", "Sire", "Dam")
                           else c("animal", "sire", "dam"))
  } else {
    .pedexport_num(ped, missing = missing)
  }

  # ---- 5. Numeric ID mapping (mirrors RENUMF90's _XrefID file) ----
  if (!use_char) {
    data.table::setattr(out, "xref", ped[, .(IndNum, Ind)])
  }

  # ---- 6. Write file if requested ----
  if (!is.null(file)) {
    .validate_export_fields(ped, software, missing, sep)
    data.table::fwrite(out, file = file, sep = sep, col.names = header,
                       quote = FALSE)
    if (!use_char) {
      xref_file <- paste0(file, ".xref")
      data.table::fwrite(attr(out, "xref"), file = xref_file, sep = sep,
                         col.names = FALSE, quote = FALSE)
      message(sprintf("Written ID mapping to: %s", xref_file))
    }
    message(sprintf("Written %d individuals to: %s", nrow(out), file))
  }

  invisible(out[])
}


# ---- Internal helpers --------------------------------------------------------

#' Ensure the ordering columns needed for export are present
#'
#' \code{tidyped(..., addnum = FALSE)} omits \code{IndNum}, \code{SireNum},
#' and \code{DamNum}, and \code{tidyped(..., addgen = FALSE)} omits
#' \code{Gen}.  These columns are optional upstream, but the export formats
#' rely on them for parent-before-offspring ordering.  When absent they are
#' reconstructed from \code{Ind}/\code{Sire}/\code{Dam} on a copy of the
#' object, leaving the caller's \code{ped} untouched.
#'
#' @param ped A complete tidyped object.
#' @return \code{ped} itself (when complete) or a copy with the missing
#'   columns filled in.
#' @keywords internal
.ensure_export_cols <- function(ped) {
  num_cols <- c("IndNum", "SireNum", "DamNum")
  has_num  <- all(num_cols %in% names(ped))
  has_gen  <- "Gen" %in% names(ped)
  if (has_num && has_gen) return(ped)

  ped <- data.table::copy(ped)
  if (!has_num) {
    # tidyped() always topologically sorts rows, so .I is already
    # parent-before-offspring; match() encodes missing parents as 0L.
    old_num_cols <- intersect(num_cols, names(ped))
    if (length(old_num_cols) > 0L) {
      ped[, (old_num_cols) := NULL]
    }
    ped[, IndNum  := .I]
    ped[, SireNum := match(Sire, Ind, nomatch = 0L)]
    ped[, DamNum  := match(Dam, Ind, nomatch = 0L)]
  }
  if (!has_gen) {
    # Without generation numbers, numeric index order is the correct
    # parent-before-offspring order for the ASReml export.
    ped[, Gen := IndNum]
  }
  ped[]
}

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

#' Build a character export table (ASReml or sommer format)
#' @param ped A complete tidyped object.
#' @param missing Character scalar for missing parents.  When \code{NA},
#'   missing parents are left as \code{NA} (the R convention used by the
#'   \code{"sommer"} format).
#' @param cols Character vector of length 3 with the output column names.
#' @return A data.table with the requested character columns.
#' @keywords internal
.pedexport_char <- function(ped, missing = "0",
                            cols = c("animal", "sire", "dam")) {
  sire_char <- ped$Sire
  dam_char  <- ped$Dam
  if (!is.na(missing)) {
    sire_char[is.na(sire_char)] <- missing
    dam_char[is.na(dam_char)]   <- missing
  }

  out <- data.table::data.table(
    animal = ped$Ind,
    sire   = sire_char,
    dam    = dam_char,
    Gen    = ped$Gen
  )
  data.table::setorderv(out, "Gen")
  out[, Gen := NULL]
  data.table::setnames(out, c("animal", "sire", "dam"), cols)
  out[]
}

#' Validate the separator supported by an export format
#' @param software Target software name.
#' @param sep Single-character field separator.
#' @return \code{NULL}, invisibly.
#' @noRd
.validate_export_sep <- function(software, sep) {
  if (software == "blupf90" && sep != " ") {
    stop("BLUPF90 pedigree files require sep = \" \".", call. = FALSE)
  }

  whitespace_formats <- c("wombat", "mtdfreml", "dmu")
  if (software %in% whitespace_formats && !(sep %in% c(" ", "\t"))) {
    stop(sprintf("%s pedigree files require a space or TAB separator.",
                 toupper(software)), call. = FALSE)
  }

  if (software == "asreml" && !(sep %in% c(" ", "\t", ","))) {
    stop("ASReml pedigree files require a space, TAB, or comma separator.",
         call. = FALSE)
  }

  invisible(NULL)
}

#' Validate unquoted text fields before writing
#' @param ped A complete tidyped object with export columns.
#' @param software Target software name.
#' @param missing Missing-parent symbol.
#' @param sep Single-character field separator.
#' @return \code{NULL}, invisibly.
#' @noRd
.validate_export_fields <- function(ped, software, missing, sep) {
  fields <- ped$Ind
  if (software %in% c("asreml", "wombat")) {
    fields <- c(fields, ped$Sire, ped$Dam, missing)
  }
  fields <- as.character(fields[!is.na(fields)])

  if (any(grepl("[[:space:]]", fields))) {
    stop("Identifiers and character missing-parent symbols must not contain ",
         "whitespace when writing an export file.", call. = FALSE)
  }

  if (!(sep %in% c(" ", "\t")) &&
      any(grepl(sep, fields, fixed = TRUE))) {
    stop("Identifiers and character missing-parent symbols must not contain ",
         "the selected separator.", call. = FALSE)
  }

  invisible(NULL)
}
