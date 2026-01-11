# Calculate inbreeding coefficients

`inbreed` function calculates the inbreeding coefficients for all
individuals in a tidied pedigree.

## Usage

``` r
inbreed(ped, ...)
```

## Arguments

- ped:

  A `tidyped` object.

- ...:

  Additional arguments passed to
  [`makeDiiF`](https://rdrr.io/pkg/nadiv/man/makeTinv.html).

## Value

A `tidyped` object with an additional column **f**.

## Details

This function takes a pedigree tidied by the
[`tidyped`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
function and calculates the inbreeding coefficients using the `makeDiiF`
function from the **nadiv** package. It prefers using numeric columns
(**IndNum**, **SireNum**, **DamNum**) if available, which is faster and
more robust.
