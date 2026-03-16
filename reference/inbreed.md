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
code based on the Meuwissen & Luo (1992) algorithm. It is the core
engine used by both `tidyped(..., inbreed = TRUE)` and
`pedmat(..., method = "f")`, ensuring consistent results across the
package. It is significantly faster than standard R implementations for
large pedigrees.

## Examples

``` r
library(visPedigree)
data(simple_ped)
ped <- tidyped(simple_ped)
ped_f <- inbreed(ped)
ped_f[f > 0, .(Ind, Sire, Dam, f)]
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#>       Ind   Sire    Dam          f
#>    <char> <char> <char>      <num>
#> 1: J5X804 J4Y326 J4E185 0.00390625
```
