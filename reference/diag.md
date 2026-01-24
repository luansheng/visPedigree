# Extract Matrix Diagonals

Extract the diagonal of a pedmat object or other matrices. This generic
function extending [`base::diag`](https://rdrr.io/r/base/diag.html) to
support `pedmat` objects.

## Usage

``` r
diag(x, ...)
```

## Arguments

- x:

  A matrix, vector, or `pedmat` object.

- ...:

  Additional arguments passed to
  [`base::diag`](https://rdrr.io/r/base/diag.html).

## Value

The diagonal vectors or a diagonal matrix.
