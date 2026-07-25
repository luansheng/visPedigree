# Efficient Pedigree Data Workflows in R

This vignette summarizes efficient day-to-day workflows for
`visPedigree` after the `tidyped` architecture updates. The goal is
simple:

1.  tidy once,
2.  reuse the resulting `tidyped` object many times,
3.  subset safely,
4.  trace candidates explicitly when pedigree completeness matters.

For basic tidying, see `tidy-pedigree`. For downstream statistics, see
`pedigree-analysis`.

## 1. Load packages and example data

``` r

library(visPedigree)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%

data(simple_ped, package = "visPedigree")
```

## 2. Tidy once, reuse many times

The most efficient workflow is to create a master `tidyped` object once
and reuse it for plotting, tracing, inbreeding, and matrix calculations.

``` r

tp_master <- tidyped(simple_ped)

class(tp_master)
#> [1] "tidyped"    "data.table" "data.frame"
is_tidyped(tp_master)
#> [1] TRUE
pedmeta(tp_master)
#> $selfing
#> [1] FALSE
#> 
#> $bisexual_parents
#> character(0)
#> 
#> $genmethod
#> [1] "top"
```

This avoids repeated validation, founder insertion, loop checking,
generation assignment, and integer re-indexing.

## 3. Fast repeated tracing from an existing `tidyped`

When the input is already a `tidyped` object and `cand` is supplied,
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
now uses a fast path. It skips the expensive global preprocessing steps
and directly traces the requested candidates.

``` r

tp_up <- tidyped(tp_master, cand = "J5X804", trace = "up", tracegen = 2)
tp_down <- tidyped(tp_master, cand = "J0Z990", trace = "down")

has_candidates(tp_up)
#> [1] TRUE
tp_up[, .(Ind, Sire, Dam, Cand)]
#> Tidy Pedigree Object
#>       Ind   Sire    Dam   Cand
#>    <char> <char> <char> <lgcl>
#> 1: J3L886   <NA>   <NA>  FALSE
#> 2: J3X697   <NA>   <NA>  FALSE
#> 3: J3Y620   <NA>   <NA>  FALSE
#> 4: J3Y771   <NA>   <NA>  FALSE
#> 5: J4E185 J3L886 J3X697  FALSE
#> 6: J4Y326 J3Y620 J3Y771  FALSE
#> 7: J5X804 J4Y326 J4E185   TRUE
```

Recommended pattern:

``` r

# expensive once
# tp_master <- tidyped(raw_ped)

# cheap many times
# tp_a <- tidyped(tp_master, cand = ids_a, trace = "up")
# tp_b <- tidyped(tp_master, cand = ids_b, trace = "all", tracegen = 3)
# tp_c <- tidyped(tp_master, cand = ids_c, trace = "down")
```

## 4. Safe `data.table` usage on `tidyped`

A `tidyped` object is also a `data.table`, so by-reference workflows
remain available.

### 4.1 Adding new columns is safe

``` r

tp_work <- copy(tp_master)
tp_work[, phenotype := seq_len(.N)]

class(tp_work)
#> [1] "tidyped"    "data.table" "data.frame"
head(tp_work[, .(Ind, phenotype)])
#>       Ind phenotype
#>    <char>     <int>
#> 1: J0C032         1
#> 2: J0C185         2
#> 3: J0C231         3
#> 4: J0C317         4
#> 5: J0C355         5
#> 6: J0C450         6
```

The `tidyped` class is preserved after `:=` operations.

### 4.2 Incomplete row subsetting now degrades safely

If row filtering removes required parents, the result is no longer a
complete pedigree. In that case the object is downgraded to a plain
`data.table` with a warning.

``` r

ped_year <- data.table(
  Ind = c("A", "B", "C", "D"),
  Sire = c(NA, NA, "A", "C"),
  Dam = c(NA, NA, "B", "B"),
  Year = c(2000, 2000, 2005, 2006)
)

tp_year <- tidyped(ped_year)
sub_dt <- tp_year[Year > 2005]
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.

class(sub_dt)
#> [1] "data.table" "data.frame"
sub_dt
#>       Ind   Sire    Dam  Year Family FamilySize   Gen    Sex IndNum SireNum
#>    <char> <char> <char> <num> <char>      <int> <int> <char>  <int>   <int>
#> 1:      D      C      B  2006    CxB          1     3   <NA>      4       3
#>    DamNum
#>     <int>
#> 1:      2
```

This behavior prevents invalid integer pedigree indices from silently
reaching C++ code.

Completeness-sensitive analyses now fail fast on such truncated subsets:

``` r

inbreed(sub_dt)
#> Error:
#> ! inbreed() requires a structurally complete pedigree. This input appears to be a row-truncated subset with missing parent records.
#> Compute on the full pedigree first, or extract a valid sub-pedigree with `tidyped(tp, cand = ids, trace = "up")`.
```

### 4.3 Use explicit tracing when you need a valid sub-pedigree

If the goal is to keep a structurally valid pedigree around focal
individuals, use candidate tracing instead of ad hoc row filtering.

``` r

valid_sub_tp <- tidyped(tp_year, cand = "D", trace = "up")

class(valid_sub_tp)
#> [1] "tidyped"    "data.table" "data.frame"
valid_sub_tp[, .(Ind, Sire, Dam, Cand)]
#> Tidy Pedigree Object
#>       Ind   Sire    Dam   Cand
#>    <char> <char> <char> <lgcl>
#> 1:      A   <NA>   <NA>  FALSE
#> 2:      B   <NA>   <NA>  FALSE
#> 3:      C      A      B  FALSE
#> 4:      D      C      B   TRUE
```

Then compute on the valid sub-pedigree and, if needed, filter the final
result back to the focal individuals:

``` r

inbreed(valid_sub_tp)[Ind == "D", .(Ind, f)]
#>       Ind     f
#>    <char> <num>
#> 1:      D  0.25
```

## 5. `splitped()` versus `pedsubpop()`

These two functions serve different purposes.

- [`splitped()`](https://luansheng.github.io/visPedigree/reference/splitped.md)
  returns the actual split pedigree objects.
- [`pedsubpop()`](https://luansheng.github.io/visPedigree/reference/pedsubpop.md)
  returns a summary table.

``` r

sub_tps <- splitped(tp_master)
length(sub_tps)
#> [1] 2
class(sub_tps[[1]])
#> [1] "tidyped"    "data.table" "data.frame"

pedsubpop(tp_master)
#>     Group     N N_Sire N_Dam N_Founder
#>    <char> <int>  <int> <int>     <int>
#> 1:    GP1    56     27    27        26
#> 2:    GP2     3      1     1         2
```

Use
[`splitped()`](https://luansheng.github.io/visPedigree/reference/splitped.md)
when you need downstream analysis on each component. Use
[`pedsubpop()`](https://luansheng.github.io/visPedigree/reference/pedsubpop.md)
when you only need the component summary.

## 6. Use accessors instead of manual attribute checks

The updated accessors are the preferred way to inspect object state.

``` r

tp_f <- inbreed(tp_master)

is_tidyped(tp_f)
#> [1] TRUE
has_inbreeding(tp_f)
#> [1] TRUE
has_candidates(tp_f)
#> [1] FALSE
pedmeta(tp_f)
#> $selfing
#> [1] FALSE
#> 
#> $bisexual_parents
#> character(0)
#> 
#> $genmethod
#> [1] "top"
```

This is preferable to hand-written checks such as `"f" %in% names(tp)`
or manual attribute access scattered throughout user code.

## 7. Recommended high-efficiency workflow

A practical pattern for large pedigrees is:

``` r

# 1. build one validated master object
# tp_master <- tidyped(raw_ped)

# 2. add analysis-specific columns in place
# tp_master[, phenotype := pheno_vector]
# tp_master[, cohort := year_vector]

# 3. extract valid candidate sub-pedigrees explicitly
# tp_sel <- tidyped(tp_master, cand = selected_ids, trace = "up", tracegen = 3)

# 4. run downstream analysis on either the full master or traced sub-pedigree
# pedstats(tp_master)
# pedmat(tp_sel)
# inbreed(tp_sel)
# visped(tp_sel)

# 5. split only when disconnected components really matter
# comps <- splitped(tp_master)
```

## 8. Practical rules of thumb

1.  Call
    [`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
    on raw pedigree data once.
2.  Reuse the resulting `tidyped` object as the master pedigree.
3.  Use `tidyped(tp_master, cand = ...)` for valid local extraction.
4.  Use ordinary row filtering only when a plain `data.table` result is
    acceptable.
5.  Use
    [`splitped()`](https://luansheng.github.io/visPedigree/reference/splitped.md)
    for actual component objects and
    [`pedsubpop()`](https://luansheng.github.io/visPedigree/reference/pedsubpop.md)
    for summaries.
6.  Use
    [`pedmeta()`](https://luansheng.github.io/visPedigree/reference/pedmeta.md),
    [`is_tidyped()`](https://luansheng.github.io/visPedigree/reference/is_tidyped.md),
    [`has_inbreeding()`](https://luansheng.github.io/visPedigree/reference/has_inbreeding.md),
    and
    [`has_candidates()`](https://luansheng.github.io/visPedigree/reference/has_candidates.md)
    to inspect object state.

These rules keep workflows fast, explicit, and structurally safe.
