# visPedigree

Tidying and Visualizing Animal Pedigrees

The `visPedigree` package provides tools to check for duplicate and
bisexual individuals, detect pedigree loops, add missing founders, sort
parents before offspring, and trace the pedigrees of specified
candidates. It generates hierarchical graphs for all individuals in a
pedigree and can handle very large datasets (\> 10,000 individuals per
generation) by compacting full-sib groups. It is particularly effective
for aquatic animal pedigrees, which often include numerous full-sib
families per generation in nucleus breeding populations.

More complex pedigree graphs and usage examples can be found in the
package vignettes.

## Installation

### From CRAN

*(Note: The package is preparing for submission and not yet available on
CRAN)*

``` r
# install.packages("visPedigree")
```

### From GitHub

You can install the development version from GitHub using the `devtools`
package:

``` r
# install.packages("devtools")
devtools::install_github("luansheng/visPedigree")
```

## Quick Start

``` r
library(visPedigree)

# Example 1: Tidy and visualize a small pedigree
# Draw the pedigree, compacting full-sib individuals and highlighting candidates.
# Note: the compacted candidates are also highlighted here.
cands <- c("Y", "Z1", "Z2")
small_ped |>
  tidyped(cand = cands) |>
  visped(compact = TRUE, highlight = cands)

# Example 2: Calculate and show inbreeding coefficients
library(data.table)
data.table(
  Ind = c("A", "B", "C", "D", "E"),
  Sire = c(NA, NA, "A", "C", "C"),
  Dam = c(NA, NA, "B", "B", "D"),
  Sex = c("male", "female", "male", "female", "male")
) |>
  tidyped(inbreed = TRUE) |>
  visped(highlight = c("E"), showf = TRUE)

# Example 3: Summarize pedigree statistics
# Get a quick overview of the pedigree structure
small_ped |>
  tidyped() |>
  summary()
```

## Citation

Luan Sheng (2026). visPedigree: Tidying and Visualizing Animal
Pedigrees. R package version 1.0.0,
<https://github.com/luansheng/visPedigree>.
