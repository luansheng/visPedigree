# Visualize Pedigree Statistics (internal)

Internal plotting backend for `plot.pedstats`. Users should call
`plot(stats, ...)` instead of this function directly.

## Usage

``` r
vispstat(x, type = c("genint", "ecg"), metric = "ECG", ...)

# S3 method for class 'pedstats'
plot(x, ...)
```

## Arguments

- x:

  A `pedstats` object returned by
  [`pedstats`](https://luansheng.github.io/visPedigree/reference/pedstats.md).

- type:

  Character. The type of plot to generate:

  - `"genint"`: Bar chart of mean generation intervals.

  - `"ecg"`: Histogram of ancestral depth metrics (ECG, FullGen, or
    MaxGen).

- metric:

  Character. Specific metric to plot when `type = "ecg"`. Options:
  `"ECG"` (default), `"FullGen"`, `"MaxGen"`.

- ...:

  Additional arguments passed to
  [`barchart`](https://rdrr.io/pkg/lattice/man/xyplot.html) or
  [`histogram`](https://rdrr.io/pkg/lattice/man/histogram.html).

## Value

A lattice plot object.

## See also

[`pedstats`](https://luansheng.github.io/visPedigree/reference/pedstats.md),
`plot.pedstats`
