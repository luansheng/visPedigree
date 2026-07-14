## Test environments
* local macOS Tahoe 26.3, R 4.5.2, Apple clang 17.0.0
* `devtools::test()`: PASS / 0 FAIL / 0 SKIP

## R CMD check results

* local `R CMD build` + `R CMD check --as-cran`
* 0 errors | 0 warnings | 0 notes (with `--no-manual`)
* `pdflatex` is not available locally, so the PDF manual check produces
  1 ERROR + 1 WARNING — these are local-environment issues only and do not
  reflect package defects. CRAN's build machines have pdflatex available.
* Tarball size is approximately 7.4 MB, above the typical 5 MB threshold.
  The size is due to (1) two vignettes (`draw-pedigree` and
  `relationship-matrix`) with pedigree visualizations and heatmaps rendered
  at 300 DPI for publication-quality output, and (2) the `big_family_size_ped`
  example dataset (473 KB) for demonstrating large-pedigree scaling.
  The compiled C++ source also contributes to the package footprint.

This is an update from CRAN version 1.8.1 to 1.9.0.

## Changes since last CRAN version (1.8.1)

This is a feature update with a restructured NEWS.md and expanded vignette
examples:

* Restructured NEWS.md with consistent section headers and concise entries
  for all historical releases.
* Added `shapeby` argument to `visped()` for choosing between sex-based
  and role-based node shape encoding.
* Added `labelvar` argument to `visped()` for custom node labels from a
  pedigree column or character vector.
* Added SVG output support in `visped()` via file names ending in `.svg`.
* Added `pedprod()` for matrix-free A x, A X, A^{-1} x, and A^{-1} X
  products directly from a pedigree, avoiding materialization of the dense
  relationship matrix.
* Expanded the `pedprod()` vignette section with group coancestry calculation
  and founder-origin decomposition examples.
* Fixed `size2` global variable declaration for R CMD check compliance.

## Downstream dependencies
None.
