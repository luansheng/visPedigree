# Internal helper to ensure ped is a tidyped object

If the object has lost its `tidyped` class (e.g., after
[`merge()`](https://rdrr.io/r/base/merge.html),
[`rbind()`](https://rdrr.io/r/base/cbind.html), or dplyr operations) but
still contains the required columns, the class is automatically restored
with an informational message. Otherwise, an error is raised guiding the
user to call
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
or
[`as_tidyped()`](https://luansheng.github.io/visPedigree/reference/as_tidyped.md).

## Usage

``` r
ensure_tidyped(ped)
```

## Arguments

- ped:

  An object expected to be a tidyped.

## Value

A valid tidyped object.
