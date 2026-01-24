# Extract Matrix Diagonals

Extract the diagonal of a pedmat object or other matrices. This generic
function delegates to [`diag`](https://rdrr.io/r/base/diag.html),
`Matrix::diag`, or specialized S3 methods.

## Usage

``` r
diag(x, ...)
```

## Arguments

- x:

  A matrix, vector, or `pedmat` object.

- ...:

  Additional arguments passed to methods.

## Value

The diagonal vectors or a diagonal matrix.
