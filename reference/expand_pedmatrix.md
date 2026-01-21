# Expand a Compact Pedigree Matrix to Full Dimensions

Restores a compact pedmatrix to its original dimensions by mapping each
individual to their family representative's values. For non-compact
matrices, returns the matrix unchanged.

## Usage

``` r
expand_pedmatrix(x)
```

## Arguments

- x:

  A pedmatrix object from
  [`pedmatrix`](https://luansheng.github.io/visPedigree/reference/pedmatrix.md).

## Value

Matrix or vector with original pedigree dimensions:

- Matrices: Row and column names set to all individual IDs

- Vectors (e.g., method="f"): Names set to all individual IDs

The result is **not** a pedmatrix object (S3 class stripped).

## Details

For compact matrices, full-siblings within the same family will have
identical relationship values in the expanded matrix because they shared
the same representative during calculation.

## See also

[`pedmatrix`](https://luansheng.github.io/visPedigree/reference/pedmatrix.md),
[`query_relationship`](https://luansheng.github.io/visPedigree/reference/query_relationship.md)

## Examples

``` r
tped <- tidyped(small_ped)

# Compact matrix
A_compact <- pedmatrix(tped, method = "A", compact = TRUE)
dim(A_compact)  # Reduced dimensions
#> [1] 23 23

# Expand to full size
A_full <- expand_pedmatrix(A_compact)
dim(A_full)  # Original dimensions restored
#> [1] 28 28

# Non-compact matrices are returned unchanged
A <- pedmatrix(tped, method = "A", compact = FALSE)
A2 <- expand_pedmatrix(A)
identical(dim(A), dim(A2))  # TRUE
#> [1] TRUE
```
