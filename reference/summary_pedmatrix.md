# Summary Statistics for Pedigree Matrices

Computes and displays summary statistics for a pedmatrix object.

## Usage

``` r
summary_pedmatrix(x)
```

## Arguments

- x:

  A pedmatrix object from
  [`pedmatrix`](https://luansheng.github.io/visPedigree/reference/pedmatrix.md).

## Value

An object of class `"summary.pedmatrix"` with statistics including
method, dimensions, compression ratio (if compact), mean relationship,
and matrix density.

## Details

Since pedmatrix objects are often S4 sparse matrices with custom
attributes, use this function instead of the generic
[`summary()`](https://rdrr.io/r/base/summary.html) to ensure proper
display of pedigree matrix statistics.

## See also

[`pedmatrix`](https://luansheng.github.io/visPedigree/reference/pedmatrix.md)

## Examples

``` r
tped <- tidyped(small_ped)
A <- pedmatrix(tped, method = "A")
summary_pedmatrix(A)
#> Summary of Pedigree Matrix (A)
#> ========================================
#> Input Size:      28  individuals
#> Calculated Size: 28  individuals
#> 
#> Matrix Properties:
#> - Mean relationship:  0.352409 
#> - Density (non-zero): 54.08%
#> ========================================
```
