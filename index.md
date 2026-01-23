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
# Draw the pedigree, compacting full-sib individuals into family nodes (squares)
# to keep the graph legible even with many offspring.
cands <- c("Y", "Z1", "Z2")
small_ped |>
  tidyped(cand = cands) |>
  visped(compact = TRUE, highlight = cands)

# Example 2: Relationship Matrices for massive pedigrees (v1.0.0+)
# Calculate A and D matrices. Use compact = TRUE to significantly speed up 
# calculations in pedigrees with large full-sib families.
mat_a <- simple_ped |> tidyped() |> pedmatrix(method = "A", compact = TRUE)
# Visualize the matrix as a heatmap
vismat(mat_a)

# Example 3: Inbreeding Coefficients
# Calculate inbreeding coefficients (f) and display them in the graph.
simple_ped |>
  tidyped(inbreed = TRUE) |>
  visped(highlight = "J5X804", showf = TRUE, compact = TRUE)

# Example 4: Pedigree Splitting and Summary
# Detect disconnected groups in the built-in simple_ped
split_list <- simple_ped |> tidyped() |> splitped()
# The result is a list of tidyped objects, one for each group
summary(split_list[[1]])
```

## Citation

Luan Sheng (2026). visPedigree: Tidying and Visualizing Animal
Pedigrees. R package version 1.0.0,
<https://github.com/luansheng/visPedigree>.
