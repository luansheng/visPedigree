# Expand a Compact Pedigree Matrix to Full Dimensions

Restores a compact pedmat to its original dimensions by mapping each
individual to their family representative's values. For non-compact
matrices, returns the matrix unchanged.

## Usage

``` r
expand_pedmat(x)
```

## Arguments

- x:

  A pedmat object from
  [`pedmat`](https://luansheng.github.io/visPedigree/reference/pedmat.md).

## Value

Matrix or vector with original pedigree dimensions:

- Matrices: Row and column names set to all individual IDs

- Vectors (e.g., method="f"): Names set to all individual IDs

The result is **not** a pedmat object (S3 class stripped).

## Details

For compact matrices, full-siblings within the same family will have
identical relationship values in the expanded matrix because they shared
the same representative during calculation.

## See also

[`pedmat`](https://luansheng.github.io/visPedigree/reference/pedmat.md),
[`query_relationship`](https://luansheng.github.io/visPedigree/reference/query_relationship.md)

## Examples

``` r
tped <- tidyped(small_ped)

# Compact matrix
A_compact <- pedmat(tped, method = "A", compact = TRUE)
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
dim(A_compact)  # Reduced dimensions
#> [1] 27 27

# Expand to full size
A_full <- expand_pedmat(A_compact)
dim(A_full)  # Original dimensions restored
#> [1] 28 28

# Non-compact matrices are returned unchanged
A <- pedmat(tped, method = "A", compact = FALSE)
A2 <- expand_pedmat(A)
identical(dim(A), dim(A2))  # TRUE
#> [1] TRUE
```
