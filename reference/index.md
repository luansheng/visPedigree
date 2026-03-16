# Package index

## Core Functions

Main tools for processing and visualizing animal pedigrees.

- [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
  : Tidy and prepare a pedigree using graph theory
- [`as_tidyped()`](https://luansheng.github.io/visPedigree/reference/as_tidyped.md)
  : Restore the tidyped class to a manipulated pedigree
- [`is_tidyped()`](https://luansheng.github.io/visPedigree/reference/is_tidyped.md)
  : Test if an object is a tidyped
- [`has_inbreeding()`](https://luansheng.github.io/visPedigree/reference/has_inbreeding.md)
  : Check whether a tidyped object contains inbreeding coefficients
- [`has_candidates()`](https://luansheng.github.io/visPedigree/reference/has_candidates.md)
  : Check whether a tidyped object contains candidate flags
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
- [`pedancestry()`](https://luansheng.github.io/visPedigree/reference/pedancestry.md)
  : Calculate Ancestry Proportions
- [`pedcontrib()`](https://luansheng.github.io/visPedigree/reference/pedcontrib.md)
  : Calculate Founder and Ancestor Contributions
- [`pedecg()`](https://luansheng.github.io/visPedigree/reference/pedecg.md)
  : Calculate Equi-Generate Coefficient
- [`pedfclass()`](https://luansheng.github.io/visPedigree/reference/pedfclass.md)
  : Summarize Inbreeding Levels
- [`pedgenint()`](https://luansheng.github.io/visPedigree/reference/pedgenint.md)
  : Calculate Generation Intervals
- [`pediv()`](https://luansheng.github.io/visPedigree/reference/pediv.md)
  : Calculate Genetic Diversity Indicators
- [`pedmat()`](https://luansheng.github.io/visPedigree/reference/pedmat.md)
  : Genetic Relationship Matrices and Inbreeding Coefficients
- [`pedmeta()`](https://luansheng.github.io/visPedigree/reference/pedmeta.md)
  : Access pedigree metadata from a tidyped object
- [`pedne()`](https://luansheng.github.io/visPedigree/reference/pedne.md)
  : Calculate Effective Population Size
- [`pedpartial()`](https://luansheng.github.io/visPedigree/reference/pedpartial.md)
  : Calculate Partial Inbreeding
- [`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md)
  : Calculate Average Additive Genetic Relationship (\\a\_{ij}\\)
- [`pedstats()`](https://luansheng.github.io/visPedigree/reference/pedstats.md)
  : Pedigree Statistics
- [`pedsubpop()`](https://luansheng.github.io/visPedigree/reference/pedsubpop.md)
  : Pedigree Subpopulations

## Datasets

Built-in pedigree examples for testing and demonstration.

- [`big_family_size_ped`](https://luansheng.github.io/visPedigree/reference/big_family_size_ped.md)
  : A large pedigree with big family sizes
- [`complex_ped`](https://luansheng.github.io/visPedigree/reference/complex_ped.md)
  : A complex pedigree
- [`deep_ped`](https://luansheng.github.io/visPedigree/reference/deep_ped.md)
  : A deep pedigree
- [`half_founder_ped`](https://luansheng.github.io/visPedigree/reference/half_founder_ped.md)
  : A pedigree with half founders
- [`inbred_ped`](https://luansheng.github.io/visPedigree/reference/inbred_ped.md)
  : A highly inbred pedigree
- [`loop_ped`](https://luansheng.github.io/visPedigree/reference/loop_ped.md)
  : A pedigree with loops
- [`simple_ped`](https://luansheng.github.io/visPedigree/reference/simple_ped.md)
  : A simple pedigree
- [`small_ped`](https://luansheng.github.io/visPedigree/reference/small_ped.md)
  : A small pedigree

## S3 Methods

Generic methods for ‘tidyped’ and ‘pedmat’ objects.

- [`` `[`( ``*`<tidyped>`*`)`](https://luansheng.github.io/visPedigree/reference/sub-.tidyped.md)
  : Subset a tidyped object
- [`plot(`*`<tidyped>`*`)`](https://luansheng.github.io/visPedigree/reference/plot.tidyped.md)
  : Plot a tidy pedigree
- [`vispstat()`](https://luansheng.github.io/visPedigree/reference/vispstat.md)
  [`plot(`*`<pedstats>`*`)`](https://luansheng.github.io/visPedigree/reference/vispstat.md)
  : Visualize Pedigree Statistics
- [`summary(`*`<tidyped>`*`)`](https://luansheng.github.io/visPedigree/reference/summary.tidyped.md)
  : Summary method for tidyped objects
- [`summary_pedmat()`](https://luansheng.github.io/visPedigree/reference/summary_pedmat.md)
  : Summary Statistics for Pedigree Matrices
- [`print(`*`<pedcontrib>`*`)`](https://luansheng.github.io/visPedigree/reference/print.pedcontrib.md)
  : Print Founder and Ancestor Contributions
- [`print(`*`<pediv>`*`)`](https://luansheng.github.io/visPedigree/reference/print.pediv.md)
  : Print Genetic Diversity Summary
- [`print(`*`<pedstats>`*`)`](https://luansheng.github.io/visPedigree/reference/print.pedstats.md)
  : Print Pedigree Statistics
- [`print(`*`<summary.tidyped>`*`)`](https://luansheng.github.io/visPedigree/reference/print.summary.tidyped.md)
  : Print method for summary.tidyped
- [`print(`*`<tidyped>`*`)`](https://luansheng.github.io/visPedigree/reference/print.tidyped.md)
  : Print method for tidyped pedigree
