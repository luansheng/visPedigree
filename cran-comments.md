## Test environments
* local macOS Sequoia 15.4, R 4.5.2
* local `devtools::test()`: 702 PASS / 0 FAIL / 0 SKIP

## R CMD check results

* local `devtools::check()`
* 0 errors | 0 warnings | 1 note
* Note: "unable to verify current time" — an environment-level timestamp
	verification issue, not a code issue.

This is an update from CRAN version 1.0.1 to 1.7.0.

## Changes since last CRAN version (1.0.1)

### Major user-visible additions
* Added population-genetic analysis workflows including:
	- `pedrel()` for average relationship trends, with both relationship and coancestry scales
	- `pedne()` for effective population size estimation
	- `pediv()` for pedigree diversity indicators (now including `GeneDiv = 1 - MeanCoan`,
		the pedigree-based retained genetic diversity)
	- `pedancestry()` and `pedpartial()` for ancestry and partial inbreeding analysis
	- `pedcontrib()` for founder and ancestor contribution analysis
	- `pedgenint()` for generation interval estimation
	- `pedecg()` for pedigree completeness / equivalent complete generations
	- `pedfclass()` for inbreeding-class summaries
	- `pedsubpop()` for pedigree subgroup summaries
	- `pedhalflife()` for information-theoretic diversity half-life analysis
* Added `tidyped` infrastructure improvements, including metadata helpers,
	safer subsetting, and class restoration tools.
* Expanded visualization workflows:
	- `vismat()` now supports direct aggregation from compact relationship matrices
		without full expansion when `by` is used
	- `vismat()` uses a representative view for pedigrees with N > 5,000 individuals,
		displaying compact representative individuals with `ID (×n)` labels to avoid
		memory overflow on very large pedigrees
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
