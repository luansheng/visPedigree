# Package index

## Core Functions

Main tools for processing and visualizing animal pedigrees.

- [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
  : Tidy and prepare a pedigree using graph theory
- [`visped()`](https://luansheng.github.io/visPedigree/reference/visped.md)
  : Visualize a tidy pedigree
- [`splitped()`](https://luansheng.github.io/visPedigree/reference/splitped.md)
  : Split Pedigree into Disconnected Groups

## Relationship Matrices

Tools for computing and analyzing genetic relationship matrices.

- [`pedmat()`](https://luansheng.github.io/visPedigree/reference/pedmat.md)
  : Genetic Relationship Matrices and Inbreeding Coefficients
- [`expand_pedmat()`](https://luansheng.github.io/visPedigree/reference/expand_pedmat.md)
  : Expand a Compact Pedigree Matrix to Full Dimensions
- [`query_relationship()`](https://luansheng.github.io/visPedigree/reference/query_relationship.md)
  : Query Relationship Coefficients from a Pedigree Matrix
- [`vismat()`](https://luansheng.github.io/visPedigree/reference/vismat.md)
  : Visualize Relationship Matrices

## Pedigree Analysis

Functions for genealogical and genetic computations.

- [`inbreed()`](https://luansheng.github.io/visPedigree/reference/inbreed.md)
  : Calculate inbreeding coefficients

## Datasets

Built-in pedigree examples for testing and demonstration.

- [`big_family_size_ped`](https://luansheng.github.io/visPedigree/reference/big_family_size_ped.md)
  : A large pedigree with big family sizes
- [`deep_ped`](https://luansheng.github.io/visPedigree/reference/deep_ped.md)
  : A deep pedigree
- [`loop_ped`](https://luansheng.github.io/visPedigree/reference/loop_ped.md)
  : A pedigree with loops
- [`simple_ped`](https://luansheng.github.io/visPedigree/reference/simple_ped.md)
  : A simple pedigree
- [`small_ped`](https://luansheng.github.io/visPedigree/reference/small_ped.md)
  : A small pedigree

## S3 Methods

Generic methods for ‘tidyped’ and ‘pedmat’ objects.

- [`plot(`*`<tidyped>`*`)`](https://luansheng.github.io/visPedigree/reference/plot.tidyped.md)
  : Plot a tidy pedigree
- [`summary(`*`<tidyped>`*`)`](https://luansheng.github.io/visPedigree/reference/summary.tidyped.md)
  : Summary method for tidyped objects
- [`summary_pedmat()`](https://luansheng.github.io/visPedigree/reference/summary_pedmat.md)
  : Summary Statistics for Pedigree Matrices
- [`print(`*`<summary.tidyped>`*`)`](https://luansheng.github.io/visPedigree/reference/print.summary.tidyped.md)
  : Print method for summary.tidyped
- [`print(`*`<tidyped>`*`)`](https://luansheng.github.io/visPedigree/reference/print.tidyped.md)
  : Print method for tidyped pedigree
