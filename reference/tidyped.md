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
  the output.

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
is TRUE. The **f** column is included if `inbreed` is TRUE.

## Details

Compared to the legacy version, this function handles cyclic pedigrees
more robustly by detecting and reporting loops, and it is generally
faster for large pedigrees due to the use of sparse graph algorithms
from the `igraph` package. Generation assignment is performed using a
topological sorting-based algorithm that ensures parents are always
placed in a generation strictly above their offspring (TGI algorithm).

## See also

[`summary.tidyped`](https://luansheng.github.io/visPedigree/reference/summary.tidyped.md)

## Examples

``` r
library(visPedigree)
library(data.table)

# Tidy a simple pedigree
tidy_ped <- tidyped(simple_ped)
head(tidy_ped)
#> Tidy Pedigree Object
#>       Ind   Sire    Dam    Sex   Gen IndNum SireNum DamNum
#>    <char> <char> <char> <char> <int>  <int>   <int>  <int>
#> 1: J0C032   <NA>   <NA> female     1      1       0      0
#> 2: J0C185   <NA>   <NA> female     1      2       0      0
#> 3: J0C231   <NA>   <NA> female     1      3       0      0
#> 4: J0C317   <NA>   <NA>   male     1      4       0      0
#> 5: J0C450   <NA>   <NA> female     1      5       0      0
#> 6: J0C561   <NA>   <NA>   male     1      6       0      0

# Trace ancestors of a specific individual within 2 generations
tidy_ped_tracegen <- tidyped(simple_ped, cand = "J5X804", trace = "up", tracegen = 2)
head(tidy_ped_tracegen)
#> Tidy Pedigree Object
#>       Ind   Sire    Dam    Sex   Gen IndNum SireNum DamNum   Cand
#>    <char> <char> <char> <char> <int>  <int>   <int>  <int> <lgcl>
#> 1: J3L886   <NA>   <NA>   male     1      1       0      0  FALSE
#> 2: J3X697   <NA>   <NA> female     1      2       0      0  FALSE
#> 3: J3Y620   <NA>   <NA>   male     1      3       0      0  FALSE
#> 4: J3Y771   <NA>   <NA> female     1      4       0      0  FALSE
#> 5: J4E185 J3L886 J3X697 female     2      5       1      2  FALSE
#> 6: J4Y326 J3Y620 J3Y771   male     2      6       3      4  FALSE

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
#> Pedigree Summary:
#> -----------------
#> Total Individuals:  81766 
#>   - Males:    185 
#>   - Females:  261 
#>   - Unknown Sex:  81320 
#> 
#> Founders (parents unknown):  447 
#> Isolated Individuals (Gen 0):  351 
#> Maximum Generation:  4 
#> Candidates Traced:  81506 
#> 
# }
```
