# Build a unified metadata list for a tidyped object

Build a unified metadata list for a tidyped object

## Usage

``` r
build_ped_meta(
  selfing = FALSE,
  bisexual_parents = character(0),
  genmethod = "top"
)
```

## Arguments

- selfing:

  Logical: whether selfing/monoecious mode was used.

- bisexual_parents:

  Character vector of IDs that appear as both sire and dam.

- genmethod:

  Character: generation assignment method ("top" or "bottom").

## Value

A named list of pedigree metadata.
