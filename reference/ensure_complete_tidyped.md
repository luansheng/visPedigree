# Internal helper to ensure ped is a complete tidyped object

Like
[`ensure_tidyped()`](https://luansheng.github.io/visPedigree/reference/ensure_tidyped.md),
but also rejects row-truncated pedigree subsets whose referenced parents
are no longer present.

## Usage

``` r
ensure_complete_tidyped(ped, fun)
```

## Arguments

- ped:

  An object expected to be a complete tidyped pedigree.

- fun:

  Character scalar. Calling function name for the error message.

## Value

A valid, structurally complete tidyped object.
