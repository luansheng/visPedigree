# Restore the tidyped class to a manipulated pedigree

Rapidly restores the `tidyped` class to a `data.table` or `data.frame`
that was previously processed by
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
but lost its class attributes due to data manipulation.

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

This helper is intended for objects that still contain the core pedigree
columns and numeric indices, but no longer inherit from `tidyped`. A
common reproducible case is
[`rbind()`](https://rdrr.io/r/base/cbind.html) on two `tidyped`
fragments, which typically returns a plain `data.table`. Converting a
`tidyped` object to a plain `data.frame` and then subsetting it also
drops the class.

Some operations, such as [`merge()`](https://rdrr.io/r/base/merge.html)
or certain dplyr workflows, may or may not preserve the `tidyped` class
depending on the versions of data.table, dplyr, and the exact method
dispatch path used in the current R session. Therefore, `as_tidyped()`
should be viewed as a safe recovery helper rather than something only
needed after one specific verb.

Typical class-loss scenarios include:

- `rbind(tped1, tped2)` — often returns plain `data.table`

- `as.data.frame(tped)[rows, ]` — returns plain `data.frame`

- manual class removal or serialization / import workflows

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

# Simulate class loss via rbind()
tp2 <- rbind(tp[1:5], tp[6:10])
class(tp2)
#> [1] "data.table" "data.frame"
# [1] "data.table" "data.frame"

# Restore the class
tp3 <- as_tidyped(tp2)
class(tp3)
#> [1] "tidyped"    "data.table" "data.frame"
# [1] "tidyped"    "data.table" "data.frame"

# It can also restore from a plain data.frame if core columns are intact
tp_df <- as.data.frame(tp)
tp4 <- tp_df[tp_df$Gen > 1, ]
class(tp4)
#> [1] "data.frame"
# [1] "data.frame"

tp5 <- as_tidyped(tp4)
class(tp5)
#> [1] "tidyped"    "data.table" "data.frame"
# [1] "tidyped"    "data.table" "data.frame"
```
