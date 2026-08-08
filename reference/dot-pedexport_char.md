# Build a character export table (ASReml or sommer format)

Build a character export table (ASReml or sommer format)

## Usage

``` r
.pedexport_char(ped, missing = "0", cols = c("animal", "sire", "dam"))
```

## Arguments

- ped:

  A complete tidyped object.

- missing:

  Character scalar for missing parents. When `NA`, missing parents are
  left as `NA` (the R convention used by the `"sommer"` format).

- cols:

  Character vector of length 3 with the output column names.

## Value

A data.table with the requested character columns.
