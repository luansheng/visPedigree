# Visualize Pedigree Statistics

`vispstat` visualizes statistics from a `pedstats` object, including
generation intervals and ancestral depth distributions.

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

  - `"genint"`: Bar chart of generation intervals (Mean ± SD).

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

[`pedstats`](https://luansheng.github.io/visPedigree/reference/pedstats.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(visPedigree)
data(simple_ped)

# Add a birth year column for generation interval calculation
simple_ped$Year <- sample(2010:2020, nrow(simple_ped), replace = TRUE)
tped <- tidyped(simple_ped)

# Calculate statistics
stats <- pedstats(tped, timevar = "Year")

# Visualize generation intervals
vispstat(stats, type = "genint")

# Visualize ancestral depth (ECG)
vispstat(stats, type = "ecg", metric = "ECG")

# Visualize fully traced generations
vispstat(stats, type = "ecg", metric = "FullGen")

# Use the plot method
plot(stats, type = "ecg", metric = "MaxGen")
} # }
```
