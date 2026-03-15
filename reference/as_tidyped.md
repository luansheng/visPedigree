# Restore the tidyped class to a manipulated pedigree

Rapidly restores the `tidyped` class to a `data.table` or `data.frame`
that was previously processed by
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
but lost its class attributes due to data manipulation (e.g.,
[`merge()`](https://rdrr.io/r/base/merge.html),
[`rbind()`](https://rdrr.io/r/base/cbind.html), or dplyr verbs).

## Usage

``` r
as_tidyped(x)
```

## Arguments

- x:

  A `data.table` or `data.frame` that was previously a tidyped object.
  It must still contain the core columns: `Ind`, `Sire`, `Dam`, `Sex`,
  `Gen`, `IndNum`, `SireNum`, `DamNum`.

## Value

A `tidyped` object.

## Details

This is a lightweight operation that only checks for the required
columns and re-attaches the class—it does **not** re-run the full
pedigree sorting, generation inference, or loop detection.

Many common R operations silently strip custom S3 class attributes:

- `merge(tped, extra)` — returns plain `data.table`

- `rbind(tped1, tped2)` — returns plain `data.table`

- `dplyr::filter(tped, ...)` — returns `tbl_df`

- `subset(tped, ...)` — returns `data.frame`

After such operations, downstream analysis functions (e.g.,
[`pedstats`](https://luansheng.github.io/visPedigree/reference/pedstats.md),
[`pedne`](https://luansheng.github.io/visPedigree/reference/pedne.md))
will either error or automatically restore the class. You can also call
`as_tidyped()` explicitly to restore the class yourself.

## See also

[`tidyped`](https://luansheng.github.io/visPedigree/reference/tidyped.md),
[`new_tidyped`](https://luansheng.github.io/visPedigree/reference/new_tidyped.md)

## Examples

``` r
library(visPedigree)
tp <- tidyped(simple_ped)
class(tp)
#> [1] "tidyped"    "data.table" "data.frame"
# [1] "tidyped"    "data.table" "data.frame"

# Simulate class loss via merge
extra <- data.table::data.table(Ind = tp$Ind[1:5], Note = "example")
tp2 <- merge(tp, extra, by = "Ind", all.x = TRUE)
class(tp2)
#> [1] "tidyped"    "data.table" "data.frame"
# [1] "data.table" "data.frame"

# Restore the class
tp3 <- as_tidyped(tp2)
class(tp3)
#> [1] "tidyped"    "data.table" "data.frame"
# [1] "tidyped"    "data.table" "data.frame"
```
