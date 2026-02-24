# visPedigree: Tidying, Analyzing, and Visualizing Animal Pedigrees

Built on graph theory and the high-performance `data.table` framework,
`visPedigree` provides a comprehensive suite of tools for processing
animal pedigrees. By modeling pedigrees as Directed Acyclic Graphs
(DAGs) using `igraph`, it ensures robust loop detection, efficient
generation assignment, and flexible ancestry tracing.

## Key Features

- **Pedigree Tidying**: Robustly handles duplicate/bisexual individuals,
  pedigree loops, and missing founders.
- **High Performance**: Optimized for massive datasets using
  `data.table` and Rcpp-based C++ implementations.
- **High-Throughput Matrix Calculation**: Calculates Additive (A),
  Dominance (D), and Additive-by-Additive (AA) relationship matrices and
  their inverses.
- **Advanced Visualization**: Generates professional vector-based
  pedigree graphs with a unique compaction algorithm for large full-sib
  families.
- **Pedigree Splitting**: Efficiently detects and splits disconnected
  sub-populations.

## Installation

### Stable version from CRAN

``` r
install.packages("visPedigree")
```

### Development version from GitHub

``` r
# install.packages("devtools")
devtools::install_github("luansheng/visPedigree", build_vignettes = TRUE)
```

## Quick Start

``` r
library(visPedigree)

# Example 1: Tidy and visualize a small pedigree
# Use compact = TRUE to condense full-sib groups into a single family node (square) 
# labeled with the group size (e.g., "2"), keeping the graph clean and legible.
cands <- c("Y", "Z1", "Z2")
small_ped |>
  tidyped(cand = cands) |>
  visped(compact = TRUE)

# Example 2: Relationship Matrices (v1.0.0+)
# Compute the additive relationship matrix (A) using high-performance C++ algorithms.
# Use compact = TRUE to accelerate calculation for pedigrees with large full-sib families.
mat_a <- simple_ped |> tidyped() |> pedmat(method = "A", compact = TRUE)
# Visualize the relationship matrix as a heatmap with histograms
vismat(mat_a)

# Example 3: Inbreeding & Highlighting
# Calculate inbreeding coefficients (f) and display them on the graph (e.g., "ID (0.05)").
# Specific individuals can be highlighted to track their inheritance paths.
simple_ped |>
  tidyped(inbreed = TRUE) |>
  visped(highlight = "J5X804", trace = "up", showf = TRUE, compact = TRUE)

# Example 4: Pedigree Splitting
# Automatically split the pedigree into independent sub-populations (connected groups).
# Each group is returned as a standalone tidyped object for separate analysis.
split_list <- simple_ped |> tidyped() |> splitped()
summary(split_list[[1]])
```

## Citation

Luan Sheng (2026). visPedigree: Tidying and Visualizing Animal
Pedigrees. R package version 1.0.1,
<https://github.com/luansheng/visPedigree>.
