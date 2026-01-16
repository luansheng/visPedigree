# Plot a tidy pedigree

Plot a tidy pedigree

## Usage

``` r
# S3 method for class 'tidyped'
plot(x, ...)
```

## Arguments

- x:

  A `tidyped` object.

- ...:

  Additional arguments passed to
  [`visped`](https://luansheng.github.io/visPedigree/reference/visped.md).

## Value

Invisibly returns a list of graph data from
[`visped`](https://luansheng.github.io/visPedigree/reference/visped.md)
(node/edge data and layout components) used to render the pedigree; the
primary result is the plot drawn on the current device.
