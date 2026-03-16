# Subset a tidyped object

Intercepts `data.table`'s `[` method for `tidyped` objects. After
subsetting, the method checks whether the result is still a valid
pedigree (all referenced parents still present). If so, `IndNum`,
`SireNum`, and `DamNum` are rebuilt and the `tidyped` class is
preserved. If the pedigree becomes structurally incomplete (missing
parent records), the result is degraded to a plain `data.table` with a
warning. Column-only selections (missing core columns) also return a
plain `data.table`.

## Usage

``` r
# S3 method for class 'tidyped'
x[...]
```

## Arguments

- x:

  A `tidyped` object.

- i, j, ...:

  Arguments passed to the `data.table` `[` method.

## Value

A `tidyped` object if the result is still a complete pedigree, otherwise
a plain `data.table`.
