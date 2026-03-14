# Pedigree Subpopulations

Summarizes pedigree subpopulations and group structure.

## Usage

``` r
pedsubpop(ped, by = NULL)
```

## Arguments

- ped:

  A `tidyped` object.

- by:

  Character. The name of the column to group by. If NULL, summarizes
  disconnected components via
  [`splitped`](https://luansheng.github.io/visPedigree/reference/splitped.md).

## Value

A `data.table` with columns:

- `Group`: Subpopulation label.

- `N`: Total individuals.

- `N_Sire`: Number of distinct sires.

- `N_Dam`: Number of distinct dams.

- `N_Founder`: Number of founders (parents unknown).

## Details

When `by = NULL`, this function is a lightweight summary wrapper around
[`splitped`](https://luansheng.github.io/visPedigree/reference/splitped.md),
returning one row per disconnected pedigree component plus an optional
`"Isolated"` row for individuals with no known parents and no offspring.
When `by` is provided, it instead summarizes the pedigree directly by
the specified column (e.g. `"Gen"`, `"Year"`, `"Breed"`).

Use `pedsubpop()` when you want a compact analytical summary table. Use
[`splitped`](https://luansheng.github.io/visPedigree/reference/splitped.md)
when you need the actual re-tidied sub-pedigree objects for downstream
plotting or analysis.

## Examples

``` r
tp <- tidyped(simple_ped)

# Summarize disconnected pedigree components
pedsubpop(tp)
#>     Group     N N_Sire N_Dam N_Founder
#>    <char> <int>  <int> <int>     <int>
#> 1:    GP1    56     27    27        26
#> 2:    GP2     3      1     1         2

# Summarize by an existing grouping variable
pedsubpop(tp, by = "Gen")
#>    Group     N N_Sire N_Dam N_Founder
#>    <int> <int>  <int> <int>     <int>
#> 1:     1    28      0     0        28
#> 2:     2    16     14    14         0
#> 3:     3     8      7     8         0
#> 4:     4     4      4     3         0
#> 5:     5     2      2     2         0
#> 6:     6     1      1     1         0
```
