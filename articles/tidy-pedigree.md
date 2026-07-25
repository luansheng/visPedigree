# Preparing and Validating Pedigree Data in R

Pedigrees are fundamental to both animal and plant breeding. They are
used to improve the accuracy of breeding value estimation, to monitor
and manage inbreeding, and to support a wide range of downstream
analyses based on ancestry and relatedness. In applied settings,
however, pedigree data are often incomplete or poorly ordered: founders
may be omitted, parents may appear after offspring, and sex annotations
may be missing or inconsistent. In addition, pedigrees are commonly
stored in a simple three-column form (individual, sire, dam), which is
convenient for storage but less convenient for checking structure,
tracing relatives, or producing readable graphical displays.

Several software tools address parts of this workflow. On Windows,
pedigraph and [pedigree
viewer](https://bkinghor.une.edu.au/pedigree.htm) provide facilities for
pedigree trimming and display. Within R, packages such as `pedigree`,
`nadiv`, and `optiSel` support pedigree preparation and analysis, while
`kinship2` can be used to draw pedigree trees. Nevertheless, large
pedigrees continue to pose practical difficulties for data cleaning,
loop detection, multi-generation tracing, and graphical display with
limited node overlap.

The [visPedigree](https://github.com/luansheng/visPedigree) package was
developed to provide an integrated workflow for pedigree tidying,
analysis, and visualization. Built on `data.table`, C++, and `igraph`,
it supports pedigree standardization, loop detection, candidate tracing,
integer pedigree construction, generation assignment, and optional
inbreeding calculation, together with scalable pedigree visualization.
The package is designed for both animal and plant pedigrees, including
selfing or monoecious mating systems. This guide focuses on the tidying
workflow and introduces the main arguments and outputs of
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md).

The main contents of this guide are as follows:

1.  [Installation of the visPedigree package](#id_1)  
2.  [The specification of pedigree format](#id_2)  
3.  [Checking and tidying pedigree](#id_3)  
    3.1 [Introduction](#id_3-1)  
    3.2 [Pedigree loop detection](#id_3-2)  
    3.3 [Tracing the pedigree of a specific individual](#id_3-3)  
    3.4 [Creating an integer pedigree](#id_3-4)  
    3.5 [Calculating inbreeding coefficients](#id_3-5)  
    3.6 [Customizing generation assignment](#id_3-6)  
    3.7 [Summarizing the pedigree](#id_3-7)  
    3.8 [Splitting large pedigrees](#id_3-8)  
    3.9 [Working with plant pedigrees and selfing](#id_3-9)

## 1. Installation of the visPedigree package

The visPedigree package can be installed from CRAN:

``` r

install.packages("visPedigree")
```

Or from GitHub:

``` r

# install.packages("devtools")
devtools::install_github("luansheng/visPedigree")
```

## 2. Pedigree format specification

The first three columns of pedigree data must be in the order of
individual, sire, and dam IDs. The column names can be customized, but
their order must remain unchanged. Individual IDs should not be coded as
`""`, `" "`, `"0"`, `*`, or `NA`; otherwise, they will be removed from
the pedigree. Missing parents should be denoted by `NA`, `0`, or `*`.
Spaces and empty strings (`""`) will also be treated as missing parents,
though this is not recommended. Additional columns, such as sex and
generation, can also be included.

## 3. Checking and tidying pedigree

### 3.1 Introduction

The pedigree can be checked and tidied through the
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
function.

This function takes a pedigree, checks for duplicates and bisexual
individuals, detects loops, adds missing founders, sorts the pedigree,
and traces candidate pedigrees.

If the `cand` parameter is provided, only those individuals and their
ancestors or descendants are retained.

Tracing direction and the number of generations can be specified using
the `trace` and `tracegen` parameters.

Virtual generations are inferred and assigned when `addgen = TRUE`.

A numeric pedigree is generated when `addnum = TRUE`.

Sex will be inferred for all individuals if sex information is missing.
If a `Sex` column is present, values should be coded as `'male'`,
`'female'`, or `NA` (unknown). Missing sex information will be inferred
from the pedigree structure where possible.

The visPedigree package comes with multiple datasets. You can check
through the following command.

``` r

data(package="visPedigree")
```

The following code displays the `simple_ped` dataset, which contains
four columns: individual, sire, dam, and sex. Missing parents are
denoted by `'NA'`, `'0'`, or `*`. Founders are not explicitly listed,
and some parents appear after their offspring in the original data.

``` r

head(simple_ped)
#>        ID   Sire    Dam    Sex
#>    <char> <char> <char> <char>
#> 1: J4Y326 J3Y620 J3Y771   male
#> 2: J1H419 J0Z938 J0Z167 female
#> 3: J2F588     NA J1Z417 female
#> 4: J1J576 J0Z938 J0Z843   male
#> 5: J1C802 J0Z333 J0C355   male
#> 6: J2Z411 J1X971 J1J134 female
tail(simple_ped)
#>        ID   Sire    Dam    Sex
#>    <char> <char> <char> <char>
#> 1: J1E852 J0Z848 J0Z624 female
#> 2: J1H604 J0C583 J0Z380 female
#> 3: J5X804 J4Y326 J4E185 female
#> 4: J1I438 J0Z990 J0Z808   male
#> 5: J2C808 J1I975 J1F266   male
#> 6: J1K462 J0C317 J0C450 female
# The number of individuals in the pedigree dataset
nrow(simple_ped)
#> [1] 31
# Individual records with missing parents
simple_ped[Sire %in% c("0", "*", "NA", NA) |
             Dam %in% c("0", "*", "NA", NA)]
#>        ID   Sire    Dam    Sex
#>    <char> <char> <char> <char>
#> 1: J2F588     NA J1Z417 female
#> 2: J1J858 J0Z060      * female
#> 3: J3X697 J2Z903      0 female
```

Example: If we incorrectly set the female `J0Z167` as the sire of
`J2F588`,
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
will detect this bisexual conflict.

``` r

x <- data.table::copy(simple_ped)
x[ID == "J2F588", Sire := "J0Z167"]
y <- tidyped(x)
#> Error:
#> ! Sex conflict detected: The following individual(s) appear as both Sire and Dam: J0Z167. This is biologically impossible. Please check and correct the pedigree data. If this is a plant pedigree with monoecious species, set selfing = TRUE.
```

The
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
function sorts the pedigree, replaces missing parents with `NA`, ensures
parents precede their offspring, and adds missing founders.

``` r

tidy_simple_ped <- tidyped(simple_ped)
head(tidy_simple_ped)
#> Tidy Pedigree Object
#>       Ind   Sire    Dam    Sex Family FamilySize   Gen IndNum SireNum DamNum
#>    <char> <char> <char> <char> <char>      <int> <int>  <int>   <int>  <int>
#> 1: J0C032   <NA>   <NA> female   <NA>          1     1      1       0      0
#> 2: J0C185   <NA>   <NA> female   <NA>          1     1      2       0      0
#> 3: J0C231   <NA>   <NA> female   <NA>          1     1      3       0      0
#> 4: J0C317   <NA>   <NA>   male   <NA>          1     1      4       0      0
#> 5: J0C355   <NA>   <NA> female   <NA>          1     1      5       0      0
#> 6: J0C450   <NA>   <NA> female   <NA>          1     1      6       0      0
tail(tidy_simple_ped)
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#>       Ind   Sire    Dam    Sex        Family FamilySize   Gen IndNum SireNum
#>    <char> <char> <char> <char>        <char>      <int> <int>  <int>   <int>
#> 1: J3X697 J2Z903   <NA> female          <NA>          1     4     54      52
#> 2: J3Y620 J2C161 J2Z411   male J2C161xJ2Z411          1     4     55      45
#> 3: J3Y771 J2G465 J2X544 female J2G465xJ2X544          1     4     56      48
#> 4: J4E185 J3L886 J3X697 female J3L886xJ3X697          1     5     57      53
#> 5: J4Y326 J3Y620 J3Y771   male J3Y620xJ3Y771          1     5     58      55
#> 6: J5X804 J4Y326 J4E185 female J4Y326xJ4E185          1     6     59      58
#>    DamNum
#>     <int>
#> 1:      0
#> 2:     51
#> 3:     49
#> 4:     54
#> 5:     56
#> 6:     57
nrow(tidy_simple_ped)
#> [1] 59
```

In the resulting `tidy_simple_ped`, founders are added with their
inferred sex, and parents are sorted before their offspring. The number
of individuals increases from 31 to 59. The columns are renamed to
`Ind`, `Sire`, and `Dam`. Missing parents are uniformly replaced with
`NA`, and
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
provides informative messages during processing. By default,
`tidy_simple_ped` includes new columns: `Gen`, `IndNum`, `SireNum`, and
`DamNum`. These can be disabled by setting `addgen = FALSE` and
`addnum = FALSE`.

If the input dataset lacks a `Sex` column, it will be automatically
added to the tidied output.

``` r

tidy_simple_ped_no_gen_num <-
  tidyped(simple_ped, addgen = FALSE, addnum = FALSE)
    head(tidy_simple_ped_no_gen_num)
#> Tidy Pedigree Object
#>       Ind   Sire    Dam    Sex Family FamilySize
#>    <char> <char> <char> <char> <char>      <int>
#> 1: J0Z938   <NA>   <NA>   male   <NA>          1
#> 2: J0Z333   <NA>   <NA>   male   <NA>          1
#> 3: J0C561   <NA>   <NA>   male   <NA>          1
#> 4: J0Z475   <NA>   <NA>   male   <NA>          1
#> 5: J0Z511   <NA>   <NA>   male   <NA>          1
#> 6: J0Z664   <NA>   <NA>   male   <NA>          1
```

Once tidied, you can use
[`data.table::fwrite()`](https://rdrr.io/pkg/data.table/man/fwrite.html)
to export the pedigree for genetic evaluation software like ASReml.

### 3.2 Pedigree loop detection

A pedigree loop occurs when an individual is its own ancestor (e.g., A
is the parent of B, B is the parent of C, and C is the parent of A).
This is a biological impossibility and a serious error in pedigree
records. The
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
function automatically detects these cycles using graph theory
algorithms. If a loop is detected, the function will stop and provide
information about the individuals involved in the loop.

The following code demonstrates what happens when a pedigree with loops
is processed:

``` r

# loop_ped contains cycles (e.g., V -> T -> R -> P -> M -> V)
# Attempting to tidy it will result in an error
try(tidyped(loop_ped))
#> Error : Pedigree error! Pedigree loops detected:
#>  M -> P -> R -> T -> V -> M
#> F -> E -> C -> A -> F
```

Detecting loops early is crucial for ensuring the integrity of genetic
evaluations.

When saving the pedigree, missing parents should typically be replaced
with `0`.

``` r

saved_ped <- data.table::copy(tidy_simple_ped)
saved_ped[is.na(Sire), Sire := "0"]
saved_ped[is.na(Dam), Dam := "0"]
data.table::fwrite(
  x = saved_ped,
  file = tempfile(fileext = ".csv"),
  sep = ",",
  quote = FALSE
)
```

### 3.3 Tracing the pedigree of a specific individual

To trace the pedigree of specific individuals, use the `cand` parameter.
This adds a `Cand` column where `TRUE` identifies the specified
candidates. If `cand` is provided, only the candidates and their
ancestors/descendants are retained.

``` r

tidy_simple_ped_J5X804_ancestors <-
  tidyped(ped = tidy_simple_ped_no_gen_num, cand = "J5X804")
  tail(tidy_simple_ped_J5X804_ancestors)
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#>       Ind   Sire    Dam    Sex        Family FamilySize IndNum SireNum DamNum
#>    <char> <char> <char> <char>        <char>      <int>  <int>   <int>  <int>
#> 1: J3X697 J2Z903   <NA> female          <NA>          1     45      43      0
#> 2: J3Y620 J2C161 J2Z411   male J2C161xJ2Z411          1     46      37     42
#> 3: J3Y771 J2G465 J2X544 female J2G465xJ2X544          1     47      40     41
#> 4: J4E185 J3L886 J3X697 female J3L886xJ3X697          1     48      44     45
#> 5: J4Y326 J3Y620 J3Y771   male J3Y620xJ3Y771          1     49      46     47
#> 6: J5X804 J4Y326 J4E185 female J4Y326xJ4E185          1     50      49     48
#>      Gen   Cand
#>    <int> <lgcl>
#> 1:     4  FALSE
#> 2:     4  FALSE
#> 3:     4  FALSE
#> 4:     5  FALSE
#> 5:     5  FALSE
#> 6:     6   TRUE
```

By default, the function traces ancestors. You can limit the number of
generations using `tracegen`. If `tracegen` is `NULL`, all available
generations are traced.

``` r

tidy_simple_ped_J5X804_ancestors_2 <-
  tidyped(ped = tidy_simple_ped_no_gen_num,
  cand = "J5X804",
  tracegen = 2)
  print(tidy_simple_ped_J5X804_ancestors_2)
#> Tidy Pedigree Object
#>       Ind   Sire    Dam    Sex        Family FamilySize IndNum SireNum DamNum
#>    <char> <char> <char> <char>        <char>      <int>  <int>   <int>  <int>
#> 1: J3L886   <NA>   <NA>   male          <NA>          1      1       0      0
#> 2: J3X697   <NA>   <NA> female          <NA>          1      2       0      0
#> 3: J3Y620   <NA>   <NA>   male          <NA>          1      3       0      0
#> 4: J3Y771   <NA>   <NA> female          <NA>          1      4       0      0
#> 5: J4E185 J3L886 J3X697 female J3L886xJ3X697          1      5       1      2
#> 6: J4Y326 J3Y620 J3Y771   male J3Y620xJ3Y771          1      6       3      4
#> 7: J5X804 J4Y326 J4E185 female J4Y326xJ4E185          1      7       6      5
#>      Gen   Cand
#>    <int> <lgcl>
#> 1:     1  FALSE
#> 2:     1  FALSE
#> 3:     1  FALSE
#> 4:     1  FALSE
#> 5:     2  FALSE
#> 6:     2  FALSE
#> 7:     3   TRUE
```

The code above traces the ancestors of `J5X804` back two generations.

To trace descendants, set `trace = 'down'`.

There are three options for the **trace** parameter:

- “up”-trace candidates’ pedigree to ancestors;
- “down”-trace candidates’ pedigree to descendants;
- “all”-trace candidaes’ pedigree to ancestors and descendants
  simultaneously.

``` r

tidy_simple_ped_J0Z990_offspring <-
  tidyped(ped = tidy_simple_ped_no_gen_num, cand = "J0Z990", trace = "down")
  print(tidy_simple_ped_J0Z990_offspring)
#> Tidy Pedigree Object
#>       Ind   Sire    Dam    Sex Family FamilySize IndNum SireNum DamNum   Gen
#>    <char> <char> <char> <char> <char>      <int>  <int>   <int>  <int> <int>
#> 1: J0Z990   <NA>   <NA>   male   <NA>          1      1       0      0     1
#> 2: J1I438 J0Z990   <NA>   male   <NA>          1      2       1      0     2
#> 3: J2G465 J1I438   <NA>   male   <NA>          1      3       2      0     3
#> 4: J3Y771 J2G465   <NA> female   <NA>          1      4       3      0     4
#> 5: J4Y326   <NA> J3Y771   male   <NA>          1      5       0      4     5
#> 6: J5X804 J4Y326   <NA> female   <NA>          1      6       5      0     6
#>      Cand
#>    <lgcl>
#> 1:   TRUE
#> 2:  FALSE
#> 3:  FALSE
#> 4:  FALSE
#> 5:  FALSE
#> 6:  FALSE
```

Tracing the descendants of `J0Z990` reveals a total of 5 individuals.

### 3.4 Creating an integer pedigree

Certain genetic evaluation programs require integer-coded pedigrees,
where individuals are numbered consecutively to facilitate the
calculation of the additive genetic relationship matrix.

By default,
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
adds `IndNum`, `SireNum`, and `DamNum` columns. This can be disabled
with `addnum = FALSE`.

``` r

tidy_simple_ped_with_int <-
  tidyped(ped = tidy_simple_ped_no_gen_num, addnum = TRUE)
head(tidy_simple_ped_with_int)
#> Tidy Pedigree Object
#>       Ind   Sire    Dam    Sex Family FamilySize   Gen IndNum SireNum DamNum
#>    <char> <char> <char> <char> <char>      <int> <int>  <int>   <int>  <int>
#> 1: J0C032   <NA>   <NA> female   <NA>          1     1      1       0      0
#> 2: J0C185   <NA>   <NA> female   <NA>          1     1      2       0      0
#> 3: J0C231   <NA>   <NA> female   <NA>          1     1      3       0      0
#> 4: J0C317   <NA>   <NA>   male   <NA>          1     1      4       0      0
#> 5: J0C355   <NA>   <NA> female   <NA>          1     1      5       0      0
#> 6: J0C450   <NA>   <NA> female   <NA>          1     1      6       0      0
```

### 3.5 Calculating inbreeding coefficients

The inbreeding coefficient (F) of each individual can be calculated
using tidyped() or inbreed() functions. There are two options to add the
inbreeding coefficients to a tidied pedigree:

1.  Set `inbreed = TRUE` in the
    [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
    function. This will calculate the inbreeding coefficients using an
    optimized C++ implementation of algorithms for large populations
    (Meuwissen & Luo, 1992; Sargolzaei et al., 2005), including the LAP
    bucket strategy described by Sargolzaei et al. (2005), and add an
    `f` column to the tidied pedigree.
2.  Or call
    [`inbreed()`](https://luansheng.github.io/visPedigree/reference/inbreed.md)
    directly on a tidied pedigree to add the `f` column.

Both options use the same high-performance engine as
`pedmat(method = "f")`, ensuring consistent results across the package.

``` r

# Create a simple inbred pedigree
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
test_ped <- data.table(
  Ind = c("A", "B", "C", "D", "E"),
  Sire = c(NA, NA, "A", "C", "C"),
  Dam = c(NA, NA, "B", "B", "D"),
  Sex = c("male", "female", "male", "female", "male")
)
# Option 1: Calculate during tidying
tidy_test <- tidyped(test_ped, inbreed = TRUE)
head(tidy_test)
#> Tidy Pedigree Object
#>       Ind   Sire    Dam    Sex Family FamilySize   Gen IndNum SireNum DamNum
#>    <char> <char> <char> <char> <char>      <int> <int>  <int>   <int>  <int>
#> 1:      A   <NA>   <NA>   male   <NA>          1     1      1       0      0
#> 2:      B   <NA>   <NA> female   <NA>          1     1      2       0      0
#> 3:      C      A      B   male    AxB          1     2      3       1      2
#> 4:      D      C      B female    CxB          1     3      4       3      2
#> 5:      E      C      D   male    CxD          1     4      5       3      4
#>        f
#>    <num>
#> 1: 0.000
#> 2: 0.000
#> 3: 0.000
#> 4: 0.250
#> 5: 0.375

# Option 2: Calculate after tidying
tidy_test <- inbreed(tidyped(test_ped))
```

### 3.6 Customizing generation assignment

Generation inference is essential for pedigree visualization.
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
provides two methods for assigning generation numbers via the
`genmethod` parameter:

- **“top” (default)**: Top-aligned (depth-based). Founders are assigned
  to Generation 1. This is the optimal scheme for most biological
  pedigrees as it ensures all founders start at the top, preventing them
  from “drifting” to lower generations if they have fewer descendants.
- **“bottom”**: Bottom-aligned (height-based). Generations are counted
  from the bottom up, aligning terminal nodes (offspring with no further
  descendants) at the highest generation number. This is useful when you
  want to show that all current populations are at the same level, or
  when introducing unrelated exogenous parents in later years.

``` r

# Default behavior (Top-Down): J2Y434 is at Gen 3
tidy_top <- tidyped(simple_ped, genmethod = "top")
tidy_top[Ind == "J2Y434"]
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#>       Ind   Sire    Dam    Sex        Family FamilySize   Gen IndNum SireNum
#>    <char> <char> <char> <char>        <char>      <int> <int>  <int>   <int>
#> 1: J2Y434 J1C802 J1H419 female J1C802xJ1H419          1     3     50      29
#>    DamNum
#>     <int>
#> 1:     34

# Bottom-Up behavior: J2Y434 is at Gen 6
tidy_bottom <- tidyped(simple_ped, genmethod = "bottom")
tidy_bottom[Ind == "J2Y434"]
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#>       Ind   Sire    Dam    Sex        Family FamilySize   Gen IndNum SireNum
#>    <char> <char> <char> <char>        <char>      <int> <int>  <int>   <int>
#> 1: J2Y434 J1C802 J1H419 female J1C802xJ1H419          1     6     58      53
#>    DamNum
#>     <int>
#> 1:     54
```

### 3.7 Summarizing the pedigree

The [`summary()`](https://rdrr.io/r/base/summary.html) method provides a
quick overview of the pedigree statistics, including the number of
individuals, sex distribution, founders, and isolated individuals. If
inbreeding coefficients have been calculated (column `f`), the summary
will also include descriptive statistics of inbreeding.

``` r

# Summarize the tidied pedigree
summary(tidy_simple_ped)
#> Pedigree Summary
#> ================
#> 
#> Total Individuals:  59 
#>   - Males:    29 (49.2%) 
#>   - Females:  30 (50.8%) 
#> 
#> Pedigree Structure:
#>   - Founders (no parents):   28 
#>   - Both parents known:      28 
#>   - Sire only known:         2 
#>   - Dam only known:          1 
#> 
#> Generation:
#>   - Maximum:  6 
#>   - Distribution:
#>       Gen 1: 28 individuals
#>       Gen 2: 16 individuals
#>       Gen 3: 8 individuals
#>       Gen 4: 4 individuals
#>       Gen 5: 2 individuals
#>       Gen 6: 1 individuals
#> 
#> Reproduction:
#>   - Individuals with offspring:  56 
#>   - Sires:  28  (Mean=1.1, Max=2 offspring)
#>   - Dams:   28  (Mean=1.0, Max=2 offspring)
#> 
#> Full-sibling Families:
#>   - Number of families:      27 
#>   - Mean family size:        1.04
#>   - Maximum family size:     2 
#>   - Top families by size:
#>       J0Z475xJ0C612: 2
#>       J0C317xJ0C450: 1
#>       J0C561xJ0C032: 1
#>       J0C583xJ0Z380: 1
#>       J0C591xJ0C231: 1
#> 
#> ================
```

### 3.8 Splitting large pedigrees

For extremely large pedigrees, it is sometimes useful to split them into
disconnected subsets or “sub-pedigrees”. The
[`splitped()`](https://luansheng.github.io/visPedigree/reference/splitped.md)
function automatically detects disconnected components (families that
share no ancestors) and splits the pedigree into a list of smaller
`tidyped` objects.

``` r

# Split the pedigree into components
sub_pedigrees <- splitped(tidy_simple_ped)

# View summary of the split result
summary(sub_pedigrees)
#> Summary of Pedigree Split
#> =========================
#> Total individuals in groups: 59 
#> Isolated individuals (Gen=0): 0 
#> Grand total: 59 
#> Number of groups:  2 
#> 
#> Size statistics:
#>   Min:     3 
#>   Max:     56 
#>   Mean:    29.5 
#>   Median:  29.5 
#> 
#> Connectivity: Pedigree contains disconnected groups

# Access a specific sub-pedigree
# first_sub <- sub_pedigrees[[1]]
```

### 3.9 Working with plant pedigrees and selfing

In plant breeding, many species are **monoecious** — a single plant
produces both male and female gametes and can serve as both the male
(sire) and female (dam) parent. Self-fertilization (selfing) is a common
practice in crops such as rice, wheat, soybean, and maize. In a
pedigree, selfing appears as the same individual listed as both Sire and
Dam for one or more offspring.

The `selfing` parameter in
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
controls how such pedigrees are handled:

- **`selfing = FALSE`** (default): Treats any individual appearing as
  both Sire and Dam as a sex conflict and returns an error. This mode is
  suitable for animal pedigrees where bisexuality is a data error.
- **`selfing = TRUE`**: Allows the same individual to be both Sire and
  Dam. These individuals are automatically assigned
  `Sex = "monoecious"`, and the pedigree metadata records this status.

The following example demonstrates a typical self-pollinating plant
breeding scenario:

``` r

library(data.table)

# A plant pedigree with self-fertilization:
# P1 (male) × P2 (female) → F1
# F1 selfs → F2_1, F2_2
# F2_1 selfs → F3
plant_ped <- data.table(
  Ind  = c("P1", "P2", "F1", "F2_1", "F2_2", "F3"),
  Sire = c(NA,   NA,   "P1", "F1",   "F1",   "F2_1"),
  Dam  = c(NA,   NA,   "P2", "F1",   "F1",   "F2_1")
)

# Without selfing = TRUE: error — F1 and F2_1 appear as both Sire and Dam
try(tidyped(plant_ped))
#> Error : Sex conflict detected: The following individual(s) appear as both Sire and Dam: F1, F2_1. This is biologically impossible. Please check and correct the pedigree data. If this is a plant pedigree with monoecious species, set selfing = TRUE.
```

Setting `selfing = TRUE` resolves the issue. The monoecious individuals
(`F1`, `F2_1`) are automatically identified, and their sex is assigned
accordingly:

``` r

tp_plant <- tidyped(plant_ped, selfing = TRUE, inbreed = TRUE)
#> Selfing mode: 2 individual(s) appear as both Sire and Dam: F1, F2_1. These will be assigned Sex = 'monoecious'.
tp_plant[]
#> Tidy Pedigree Object
#> Indices: <Sire>, <Dam>
#>       Ind   Sire    Dam    Family FamilySize   Gen        Sex IndNum SireNum
#>    <char> <char> <char>    <char>      <int> <int>     <char>  <int>   <int>
#> 1:     P1   <NA>   <NA>      <NA>          1     1       male      1       0
#> 2:     P2   <NA>   <NA>      <NA>          1     1     female      2       0
#> 3:     F1     P1     P2     P1xP2          1     2 monoecious      3       1
#> 4:   F2_1     F1     F1     F1xF1          2     3 monoecious      4       3
#> 5:   F2_2     F1     F1     F1xF1          2     3       <NA>      5       3
#> 6:     F3   F2_1   F2_1 F2_1xF2_1          1     4       <NA>      6       4
#>    DamNum     f
#>     <int> <num>
#> 1:      0  0.00
#> 2:      0  0.00
#> 3:      2  0.00
#> 4:      3  0.50
#> 5:      3  0.50
#> 6:      4  0.75

# Check metadata
attr(tp_plant, "ped_meta")
#> $selfing
#> [1] TRUE
#> 
#> $bisexual_parents
#> [1] "F1"   "F2_1"
#> 
#> $genmethod
#> [1] "top"
```

The [`summary()`](https://rdrr.io/r/base/summary.html) method reports
the count of monoecious individuals:

``` r

summary(tp_plant)
#> Pedigree Summary
#> ================
#> 
#> Total Individuals:  6 
#>   - Males:    1 (16.7%) 
#>   - Females:  1 (16.7%) 
#>   - Monoecious:  2 (33.3%) 
#>   - Unknown:  2 (33.3%) 
#> 
#> Pedigree Structure:
#>   - Founders (no parents):   2 
#>   - Both parents known:      4 
#> 
#> Generation:
#>   - Maximum:  4 
#>   - Distribution:
#>       Gen 1: 2 individuals
#>       Gen 2: 1 individuals
#>       Gen 3: 2 individuals
#>       Gen 4: 1 individuals
#> 
#> Reproduction:
#>   - Individuals with offspring:  4 
#>   - Sires:  3  (Mean=1.3, Max=2 offspring)
#>   - Dams:   3  (Mean=1.3, Max=2 offspring)
#> 
#> Full-sibling Families:
#>   - Number of families:      3 
#>   - Mean family size:        1.33
#>   - Maximum family size:     2 
#>   - Top families by size:
#>       F1xF1: 2
#>       F2_1xF2_1: 1
#>       P1xP2: 1
#> 
#> Inbreeding Coefficients:
#>   - All individuals:
#>       Mean = 0.2917, Min = 0.0000, Max = 0.7500
#>       Inbred (f > 0): 3 (50.0%)
#> 
#> ================
```

In the tidied output: - **`P1`** → `male` (only appears as Sire) -
**`P2`** → `female` (only appears as Dam) - **`F1`** → `monoecious`
(appears as both Sire and Dam for `F2_1`, `F2_2`) - **`F2_1`** →
`monoecious` (appears as both Sire and Dam for `F3`)

Self-fertilization produces inbred offspring: the inbreeding coefficient
*f* = 0.5 for individuals derived from a single generation of selfing,
and increases with successive selfing generations (*f* = 0.75 for `F3`).

``` r

# Inbreeding coefficients
tp_plant[, .(Ind, Sex, Gen, f)]
#>       Ind        Sex   Gen     f
#>    <char>     <char> <int> <num>
#> 1:     P1       male     1  0.00
#> 2:     P2     female     1  0.00
#> 3:     F1 monoecious     2  0.00
#> 4:   F2_1 monoecious     3  0.50
#> 5:   F2_2       <NA>     3  0.50
#> 6:     F3       <NA>     4  0.75
```

See
[`vignette("draw-pedigree", package = "visPedigree")`](https://luansheng.github.io/visPedigree/articles/draw-pedigree.md)
for visualizing plant pedigrees — monoecious individuals are drawn as
**hexagons** with teal edges and fill.

**Tip:** When working with large plant breeding datasets where selfing
is expected, always set `selfing = TRUE` to avoid false sex-conflict
errors. The `selfing` attribute is preserved through
[`splitped()`](https://luansheng.github.io/visPedigree/reference/splitped.md)
and other downstream operations, ensuring consistent handling across the
entire workflow.

------------------------------------------------------------------------

**See Also:** -
[`vignette("draw-pedigree", package = "visPedigree")`](https://luansheng.github.io/visPedigree/articles/draw-pedigree.md) -
[`vignette("relationship-matrix", package = "visPedigree")`](https://luansheng.github.io/visPedigree/articles/relationship-matrix.md)

## References

- Meuwissen, T. H. E., & Luo, Z. (1992). Computing inbreeding
  coefficients in large populations. *Genetics Selection Evolution*, 24,
  305-313.
- Sargolzaei, M., Iwaisaki, H., & Colleau, J.-J. (2005). A fast
  algorithm for computing inbreeding coefficients in large populations.
  *Journal of Animal Breeding and Genetics*, 122, 325-331.
