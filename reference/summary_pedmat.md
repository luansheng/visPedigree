# Summary Statistics for Pedigree Matrices

Computes and displays summary statistics for a pedmat object.

## Usage

``` r
summary_pedmat(x)
```

## Arguments

- x:

  A pedmat object from
  [`pedmat`](https://luansheng.github.io/visPedigree/reference/pedmat.md).

## Value

An object of class `"summary.pedmat"` with statistics including method,
dimensions, compression ratio (if compact), mean relationship, and
matrix density.

## Details

Since pedmat objects are often S4 sparse matrices with custom
attributes, use this function instead of the generic
[`summary()`](https://rdrr.io/r/base/summary.html) to ensure proper
display of pedigree matrix statistics.

## See also

[`pedmat`](https://luansheng.github.io/visPedigree/reference/pedmat.md)

## Examples

``` r
tped <- tidyped(small_ped)
A <- pedmat(tped, method = "A")
summary_pedmat(A)
#> Summary of Pedigree Matrix (A)
#> ========================================
#> Input Size:      28  individuals
#> Calculated Size: 28  individuals
#> 
#> Matrix Properties:
#> - Mean off-diagonal relationship:  0.136088 
#> - Density (non-zero): 54.08%
#> ========================================
```
