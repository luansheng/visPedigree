## Test environments
* local macOS Tahoe 26.6.2, R 4.5.2 (2025-10-31), aarch64-apple-darwin20,
  Apple clang 17.0.0 (clang-1700.0.13.5)
* `R CMD check --as-cran --no-manual`: 0 errors | 0 warnings | 0 notes
* testthat: FAIL 0 | WARN 2 | SKIP 1 | PASS 1471

## R CMD check results

Ran locally with `R CMD build` + `R CMD check --as-cran`.

With `--no-manual` the check is clean: 0 errors | 0 warnings | 0 notes.

Without `--no-manual` there are 1 ERROR, 1 WARNING and 2 NOTEs, all of them
artifacts of this local machine rather than package defects:

* `pdflatex` is not installed locally, so the PDF version of the manual cannot
  be built (1 ERROR + 1 WARNING) and the failed run leaves
  `visPedigree-manual.tex` in the check directory (1 NOTE). CRAN's build
  machines have pdflatex.
* The HTML manual check is skipped because the local `tidy` is older than
  required and the `V8` package is unavailable (1 NOTE).

The two test warnings are intentional: the package warns when subsetting a
`tidyped` object removes parent records, and two tests exercise that path
without wrapping it in `expect_warning()`.

This is an update from CRAN version 1.9.0 to 1.10.1.

## Changes since last CRAN version (1.9.0)

This release adds `pedexport()` for writing pedigrees to breeding-software
formats, plus follow-up fixes to the new formats.

### New features

* `pedexport()` converts a `tidyped` pedigree into the input format of common
  animal and plant breeding programs: BLUPF90, ASReml, Echidna, WOMBAT,
  MTDFREML, DMU, a generic numeric layout, and an in-memory `sommer` table.
  Character formats keep the original IDs; numeric formats renumber the
  pedigree and carry an `xref` mapping back to the original IDs, both as an
  attribute and as a `<file>.xref` file, mirroring RENUMF90's `_XrefID`.
* `pedexport(software = "hiblup")` writes HIBLUP `--pedigree` files.
* Export validation rejects malformed separators, invalid file paths, and
  unquoted identifiers or missing-parent symbols that would split into extra
  fields.

### Bug fixes

* `pedexport(software = "wombat")` now uses the integer layout required by the
  WOMBAT manual (section 6.3): three integer columns, offspring codes
  numerically larger than either parent, and unknown parents coded `0`. This
  replaces the 1.10.0 character-ID layout, which the WOMBAT program itself
  rejected.

### Documentation

* Documented the `sommer` workflow using
  `pedmat(..., method = "A", sparse = FALSE)`.

## Downstream dependencies
None.
