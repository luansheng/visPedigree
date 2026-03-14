# Tidy and prepare a pedigree using graph theory

This function takes a pedigree, checks for duplicate and bisexual
individuals, detects pedigree loops using graph theory, adds missing
founders, assigns generation numbers, sorts the pedigree, and traces the
pedigree of specified candidates. If the `cand` parameter contains
individual IDs, only those individuals and their ancestors or
descendants are retained. Tracing direction and the number of
generations can be specified using the `trace` and `tracegen`
parameters.

## Usage

``` r
tidyped(
  ped,
  cand = NULL,
  trace = "up",
  tracegen = NULL,
  addgen = TRUE,
  addnum = TRUE,
  inbreed = FALSE,
  selfing = FALSE,
  genmethod = "top",
  ...
)
```

## Arguments

- ped:

  A data.table or data frame containing the pedigree. The first three
  columns must be **individual**, **sire**, and **dam** IDs. Additional
  columns, such as sex or generation, can be included. Column names can
  be customized, but their order must remain unchanged. Individual IDs
  should not be coded as "", " ", "0", "\*", or "NA"; otherwise, they
  will be removed. Missing parents should be denoted by "NA", "0", or
  "\*". Spaces and empty strings ("") are also treated as missing
  parents but are not recommended.

- cand:

  A character vector of individual IDs, or NULL. If provided, only the
  candidates and their ancestors/descendants are retained.

- trace:

  A character value specifying the tracing direction: "**up**",
  "**down**", or "**all**". "up" traces ancestors; "down" traces
  descendants; "all" traces the union of ancestors and descendants. This
  parameter is only used if `cand` is not NULL. Default is "up".

- tracegen:

  An integer specifying the number of generations to trace. This
  parameter is only used if `trace` is not NULL. If NULL or 0, all
  available generations are traced.

- addgen:

  A logical value indicating whether to generate generation numbers.
  Default is TRUE, which adds a **Gen** column to the output.

- addnum:

  A logical value indicating whether to generate a numeric pedigree.
  Default is TRUE, which adds **IndNum**, **SireNum**, and **DamNum**
  columns to the output.

- inbreed:

  A logical value indicating whether to calculate inbreeding
  coefficients. Default is FALSE. If TRUE, an **f** column is added to
  the output. This uses the same optimized engine as
  `pedmat(..., method = "f")`.

- selfing:

  A logical value indicating whether to allow the same individual to
  appear as both sire and dam. This is common in plant breeding
  (monoecious species) where the same plant can serve as both male and
  female parent. If TRUE, individuals appearing in both the Sire and Dam
  columns will be assigned Sex = "monoecious" instead of triggering an
  error. Default is FALSE.

- genmethod:

  A character value specifying the generation assignment method:
  "**top**" or "**bottom**". "top" (top-aligned) assigns generations
  from parents to offspring, starting founders at Gen 1. "bottom"
  (bottom-aligned) assigns generations from offspring to parents,
  aligning terminal nodes at the bottom. Default is "top".

- ...:

  Additional arguments passed to
  [`inbreed`](https://luansheng.github.io/visPedigree/reference/inbreed.md).

## Value

A `tidyped` object (which inherits from `data.table`). Individual, sire,
and dam ID columns are renamed to **Ind**, **Sire**, and **Dam**.
Missing parents are replaced with **NA**. The **Sex** column contains
"male", "female", or NA. The **Cand** column is included if `cand` is
not NULL. The **Gen** column is included if `addgen` is TRUE. The
**IndNum**, **SireNum**, and **DamNum** columns are included if `addnum`
is TRUE. The **Family** and **FamilySize** columns identify full-sibling
families (e.g., "A x B" for offspring of sire A and dam B). The **f**
column is included if `inbreed` is TRUE.

## Details

Compared to the legacy version, this function handles cyclic pedigrees
more robustly by detecting and reporting loops, and it is generally
faster for large pedigrees due to the use of sparse graph algorithms
from the `igraph` package. Generation assignment can be performed using
either a top-down approach (default, aligning founders at the top) or a
bottom-up approach (aligning terminal nodes at the bottom).

## See also

[`summary.tidyped`](https://luansheng.github.io/visPedigree/reference/summary.tidyped.md)
for summarizing tidyped objects
[`visped`](https://luansheng.github.io/visPedigree/reference/visped.md)
for visualizing pedigree structure
[`pedmat`](https://luansheng.github.io/visPedigree/reference/pedmat.md)
for computing relationship matrices
[`vismat`](https://luansheng.github.io/visPedigree/reference/vismat.md)
for visualizing relationship matrices
[`splitped`](https://luansheng.github.io/visPedigree/reference/splitped.md)
for splitting pedigree into connected components
[`inbreed`](https://luansheng.github.io/visPedigree/reference/inbreed.md)
for calculating inbreeding coefficients

## Examples

``` r
library(visPedigree)
library(data.table)

# Tidy a simple pedigree
tidy_ped <- tidyped(simple_ped)
head(tidy_ped)
#> Tidy Pedigree Object
#>       Ind   Sire    Dam    Sex Family FamilySize   Gen IndNum SireNum DamNum
#>    <char> <char> <char> <char> <char>      <int> <int>  <int>   <int>  <int>
#> 1: J0C032   <NA>   <NA> female   <NA>          1     1      1       0      0
#> 2: J0C185   <NA>   <NA> female   <NA>          1     1      2       0      0
#> 3: J0C231   <NA>   <NA> female   <NA>          1     1      3       0      0
#> 4: J0C317   <NA>   <NA>   male   <NA>          1     1      4       0      0
#> 5: J0C355   <NA>   <NA> female   <NA>          1     1      5       0      0
#> 6: J0C450   <NA>   <NA> female   <NA>          1     1      6       0      0

# Trace ancestors of a specific individual within 2 generations
tidy_ped_tracegen <- tidyped(simple_ped, cand = "J5X804", trace = "up", tracegen = 2)
head(tidy_ped_tracegen)
#> Tidy Pedigree Object
#>       Ind   Sire    Dam    Sex        Family FamilySize   Gen IndNum SireNum
#>    <char> <char> <char> <char>        <char>      <int> <int>  <int>   <int>
#> 1: J3L886   <NA>   <NA>   male          <NA>          1     1      1       0
#> 2: J3X697   <NA>   <NA> female          <NA>          1     1      2       0
#> 3: J3Y620   <NA>   <NA>   male          <NA>          1     1      3       0
#> 4: J3Y771   <NA>   <NA> female          <NA>          1     1      4       0
#> 5: J4E185 J3L886 J3X697 female J3L886xJ3X697          1     2      5       1
#> 6: J4Y326 J3Y620 J3Y771   male J3Y620xJ3Y771          1     2      6       3
#>    DamNum   Cand
#>     <int> <lgcl>
#> 1:      0  FALSE
#> 2:      0  FALSE
#> 3:      0  FALSE
#> 4:      0  FALSE
#> 5:      2  FALSE
#> 6:      4  FALSE

# Trace both ancestors and descendants for multiple candidates
# This is highly optimized and works quickly even on 100k+ individuals
cand_list <- c("J5X804", "J3Y620")
tidy_ped_all <- tidyped(simple_ped, cand = cand_list, trace = "all")

# Check for loops (will error if loops exist)
try(tidyped(loop_ped))
#> Error : Pedigree error! Pedigree loops detected:
#>  M -> P -> R -> T -> V -> M
#> F -> E -> C -> A -> F

# Example with a large pedigree: extract 2 generations of ancestors for 2007 candidates
cand_2007 <- big_family_size_ped[Year == 2007, Ind]
# \donttest{
tidy_big <- tidyped(big_family_size_ped, cand = cand_2007, trace = "up", tracegen = 2)
summary(tidy_big)
#> Pedigree Summary
#> ================
#> 
#> Total Individuals:  81766 
#>   - Males:    185 (0.2%) 
#>   - Females:  261 (0.3%) 
#>   - Unknown:  81320 (99.5%) 
#> 
#> Pedigree Structure:
#>   - Founders (no parents):   447 
#>   - Both parents known:      81319 
#>   - Isolated (Gen 0):        351 
#> 
#> Generation:
#>   - Maximum:  4 
#>   - Distribution:
#>       Gen 0: 351 individuals
#>       Gen 1: 96 individuals
#>       Gen 2: 164 individuals
#>       Gen 3: 44146 individuals
#>       Gen 4: 37009 individuals
#> 
#> Reproduction:
#>   - Individuals with offspring:  446 
#>   - Sires:  185  (Mean=439.6, Max=930 offspring)
#>   - Dams:   261  (Mean=311.6, Max=469 offspring)
#> 
#> Full-sibling Families:
#>   - Number of families:      261 
#>   - Mean family size:        311.57
#>   - Maximum family size:     469 
#>   - Top families by size:
#>       6040x6Z30: 469
#>       6040x6Z3Z: 461
#>       6007x6074: 459
#>       60Z6x6089: 459
#>       60Z8x609Y: 455
#> 
#> Candidates Traced:  81506 
#> 
#> ================
# }
```
