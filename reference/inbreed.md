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
function and calculates the inbreeding coefficients using an optimized
C++ implementation of the Sargolzaei & Iwaisaki (2005) LAP (Longest
Ancestral Path) bucket algorithm. This method is the fastest known
direct algorithm for computing all inbreeding coefficients: it replaces
the O(\\N^2\\) linear scan of Meuwissen & Luo (1992) with O(1) bucket
pops and selective ancestor clearing, giving \\O(\sum m_i)\\ total work
where \\m_i\\ is the number of distinct ancestors of individual \\i\\.
At \\N = 1{,}000{,}000\\, the kernel completes in approximately 0.12 s —
over 10\\\times\\ faster than the previous Meuwissen & Luo (1992)
implementation and on par with the pedigreemm reference C implementation
of the same algorithm. It is the core engine used by both
`tidyped(..., inbreed = TRUE)` and `pedmat(..., method = "f")`, ensuring
consistent results across the package.

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
