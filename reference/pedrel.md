# Calculate Average Additive Genetic Relationship (\\a\_{ij}\\)

Computes the average pairwise additive genetic relationship coefficients
(\\a\_{ij}\\) within cohorts or groups. The relationship \\a\_{ij}\\ is
defined as twice the coancestry coefficient (\\f\_{ij}\\), representing
the expected proportion of genes shared by descent (e.g., 0.5 for full
siblings).

## Usage

``` r
pedrel(ped, by = "Gen", reference = NULL, compact = FALSE)
```

## Arguments

- ped:

  A `tidyped` object.

- by:

  Character. The column name to group by (e.g., "Year", "Breed",
  "Generation").

- reference:

  Character vector. An optional vector of reference individual IDs to
  calculate relationships for. If provided, only individuals matching
  these IDs in each group will be used. Default is NULL (use all
  individuals in the group).

- compact:

  Logical. Whether to use compact representation for large families to
  save memory. Recommended when pedigree size exceeds 25,000. Default is
  FALSE.

## Value

A `data.table` with columns:

- A grouping identifier column, named after the `by` parameter (e.g.,
  `Gen`, `Year`).

- `NTotal`: Total number of individuals in the group.

- `NUsed`: Number of individuals used in calculation (could be subset by
  reference).

- `MeanRel`: Average of off-diagonal elements in the Additive
  Relationship (A) matrix for this group (\\a\_{ij} = 2f\_{ij}\\).
  Returns NA if the group has fewer than 2 individuals.

## Examples

``` r
# \donttest{
library(data.table)
# Use the sample dataset and simulate a birth year
tp <- tidyped(small_ped)
tp$Year <- 2010 + tp$Gen

# Example 1: Calculate average relationship grouped by Generation (default)
rel_by_gen <- pedrel(tp, by = "Gen")
print(rel_by_gen)
#>      Gen NTotal NUsed   MeanRel
#>    <int>  <int> <int>     <num>
#> 1:     1      9     9 0.0000000
#> 2:     2      5     5 0.2000000
#> 3:     3      7     7 0.1547619
#> 4:     4      3     3 0.1145833
#> 5:     5      2     2 0.0468750
#> 6:     6      2     2 0.5507812

# Example 2: Calculate average relationship grouped by Year
rel_by_year <- pedrel(tp, by = "Year")
print(rel_by_year)
#>     Year NTotal NUsed   MeanRel
#>    <num>  <int> <int>     <num>
#> 1:  2011      9     9 0.0000000
#> 2:  2012      5     5 0.2000000
#> 3:  2013      7     7 0.1547619
#> 4:  2014      3     3 0.1145833
#> 5:  2015      2     2 0.0468750
#> 6:  2016      2     2 0.5507812

# Example 3: Filter calculations with a reference list in a chosen group
candidates <- c("N", "O", "P", "Q", "T", "U", "V", "X", "Y")
rel_subset <- pedrel(tp, by = "Gen", reference = candidates)
#> Warning: Group '3' has less than 2 individuals after applying 'reference', returning NA_real_.
#> Warning: Group '6' has less than 2 individuals after applying 'reference', returning NA_real_.
print(rel_subset)
#>      Gen NTotal NUsed  MeanRel
#>    <int>  <int> <int>    <num>
#> 1:     1      9     2 0.000000
#> 2:     2      5     2 0.500000
#> 3:     3      7     1       NA
#> 4:     4      3     2 0.156250
#> 5:     5      2     2 0.046875
#> 6:     6      2     0       NA
# }
```
