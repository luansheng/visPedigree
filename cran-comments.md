## Test environments
* local macOS Sequoia 15.4, R 4.5.2
* local `devtools::test()`: 698 PASS / 0 FAIL / 0 SKIP

## R CMD check results

* local `devtools::check(args = "--no-tests")`
* 0 errors | 0 warnings | 0 notes

This is an update from CRAN version 1.0.1 to 1.6.2.

## Changes since last CRAN version (1.0.1)

### Major user-visible additions
* Added population-genetic analysis workflows including:
	- `pedrel()` for average relationship trends, with both relationship and coancestry scales
	- `pedne()` for effective population size estimation
	- `pediv()` for pedigree diversity indicators
	- `pedancestry()` and `pedpartial()` for ancestry and partial inbreeding analysis
	- `pedcontrib()` for founder and ancestor contribution analysis
	- `pedgenint()` for generation interval estimation
	- `pedecg()` for pedigree completeness / equivalent complete generations
	- `pedfclass()` for inbreeding-class summaries
	- `pedsubpop()` for pedigree subgroup summaries
	- `pedhalflife()` for information-theoretic diversity half-life analysis
* Added robust `tidyped` infrastructure improvements, including metadata helpers,
	safer subsetting, and class restoration tools.
* Expanded visualization workflows:
	- `vismat()` now supports direct aggregation from compact relationship matrices
		without full expansion when `by` is used
	- large compact matrices can now be shown in representative view when full
		expansion is not practical
	- `plot.pedstats()` is now the main user-facing visualization route for pedigree
		statistics
* Added/expanded package vignettes covering tidy pedigree workflows, relationship
	matrices, pedigree analysis, and drawing pedigrees.

### Notable bug fixes and quality improvements
* Strengthened fail-fast behavior for incomplete pedigrees in completeness-sensitive
	analyses.
* Fixed false-positive internal warnings caused by subsetting `tidyped` objects in
	downstream analysis helpers.
* Improved compact matrix handling, including sibling off-diagonal correction and
	representative-view visualization for large pedigrees.
* Standardized and expanded tests, examples, and documentation across the package.

## Downstream dependencies
None.
