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

  Additional arguments (currently ignored).

## Value

A `tidyped` object with an additional column **f**.

## Details

This function takes a pedigree tidied by the
[`tidyped`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
function and calculates the inbreeding coefficients using optimized C++
code based on the Meuwissen & Luo (1992) algorithm. It is significantly
faster than standard R implementations for large pedigrees.
