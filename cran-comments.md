## Test environments
* local macOS Tahoe 26.3, R 4.5.2
* local `devtools::test()`: 715 PASS / 0 FAIL / 0 SKIP

## R CMD check results

* local `devtools::check(cran = TRUE, manual = TRUE, error_on = "never")`
* 0 errors | 0 warnings | 0 notes
* `pdflatex` is not available locally, so the check ran as `--no-manual --as-cran`.

This is an update from CRAN version 1.0.1 to 1.8.1.

## Changes since last CRAN version (1.0.1)

This is a substantial feature update. Main changes include:

* Added new pedigree-analysis workflows, including `pedrel()`, `pedne()`,
	`pediv()`, `pedancestry()`, `pedpartial()`, `pedcontrib()`, `pedgenint()`,
	`pedecg()`, `pedfclass()`, `pedsubpop()`, and `pedhalflife()`.
* Extended the `tidyped` infrastructure with stricter structural validation,
	safer subsetting behavior, and class/metadata restoration helpers.
* Improved visualization workflows, including compact matrix aggregation in
	`vismat()`, representative large-pedigree views, and `plot.pedstats()` as the
	main user-facing route for pedigree statistics plots.
* Added and expanded vignettes for tidy pedigree workflows, relationship
	matrices, pedigree analysis, and pedigree drawing.
* Fixed several correctness and robustness issues, especially around incomplete
	pedigrees, internal `tidyped` subsetting, and compact matrix handling.

## Downstream dependencies
None.
