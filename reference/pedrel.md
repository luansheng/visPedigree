# Calculate Mean Relationship or Coancestry Within Groups

Computes either the average pairwise additive genetic relationship
coefficients (\\a\_{ij}\\) within cohorts, or the corrected population
mean coancestry used for pedigree-based diversity summaries.

## Usage

``` r
pedrel(
  ped,
  by = "Gen",
  reference = NULL,
  compact = FALSE,
  scale = c("relationship", "coancestry")
)
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

  Logical. Retained for backward compatibility. It is ignored because
  `pedrel()` now uses matrix-free pedigree products for all
  calculations.

- scale:

  Character. One of `"relationship"` or `"coancestry"`. `"relationship"`
  returns the pairwise off-diagonal mean additive relationship (current
  `pedrel()` behavior). `"coancestry"` returns the corrected population
  mean coancestry used for pedigree-based diversity calculations.

## Value

A `data.table` with columns:

- A grouping identifier column, named after the `by` parameter (e.g.,
  `Gen`, `Year`).

- `NTotal`: Total number of individuals in the group.

- `NUsed`: Number of individuals used in calculation (could be subset by
  reference).

- `MeanRel`: Present when `scale = "relationship"`; average of
  off-diagonal elements in the Additive Relationship (A) matrix for this
  group (\\a\_{ij} = 2f\_{ij}\\).

- `MeanCoan`: Present when `scale = "coancestry"`; diagonal-corrected
  population mean coancestry for this group.

- `Status`: `"ok"`, `"skipped"`, or `"failed"`.

- `Message`: Empty for successful groups; otherwise a diagnostic
  explaining why `NA` was returned.

## Details

Let \\g\\ be a zero-one indicator vector for the \\N\\ selected
individuals and let \\S = g^\prime A g\\. When `scale = "relationship"`,
the returned value is \$\$\bar{a}\_{off} = \frac{S - \sum\_{i:g_i=1}(1 +
F_i)}{N(N - 1)},\$\$ the mean of the off-diagonal additive relationship
coefficients. When `scale = "coancestry"`, the returned value is
\$\$\bar{C} = \frac{S}{2N^2},\$\$ which includes both pairwise
coancestry and diagonal self-coancestry. This is equivalent to
\$\$\bar{C} = \frac{N - 1}{N} \cdot \frac{\bar{a}\_{off}}{2} + \frac{1 +
\bar{F}}{2N}.\$\$ The implementation computes \\A g\\ directly from the
pedigree and does not construct the dense relationship matrix \\A\\.

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
#>      Gen NTotal NUsed   MeanRel Status Message
#>    <int>  <int> <int>     <num> <char>  <char>
#> 1:     1      9     9 0.0000000     ok        
#> 2:     2      5     5 0.2000000     ok        
#> 3:     3      7     7 0.1547619     ok        
#> 4:     4      3     3 0.1145833     ok        
#> 5:     5      2     2 0.0468750     ok        
#> 6:     6      2     2 0.5507812     ok        

# Example 2: Calculate average relationship grouped by Year
rel_by_year <- pedrel(tp, by = "Year")
print(rel_by_year)
#>     Year NTotal NUsed   MeanRel Status Message
#>    <num>  <int> <int>     <num> <char>  <char>
#> 1:  2011      9     9 0.0000000     ok        
#> 2:  2012      5     5 0.2000000     ok        
#> 3:  2013      7     7 0.1547619     ok        
#> 4:  2014      3     3 0.1145833     ok        
#> 5:  2015      2     2 0.0468750     ok        
#> 6:  2016      2     2 0.5507812     ok        

# Example 3: Calculate corrected mean coancestry
coan_by_gen <- pedrel(tp, by = "Gen", scale = "coancestry")
print(coan_by_gen)
#>      Gen NTotal NUsed   MeanCoan Status Message
#>    <int>  <int> <int>      <num> <char>  <char>
#> 1:     1      9     9 0.05555556     ok        
#> 2:     2      5     5 0.18000000     ok        
#> 3:     3      7     7 0.13775510     ok        
#> 4:     4      3     3 0.20833333     ok        
#> 5:     5      2     2 0.27148438     ok        
#> 6:     6      2     2 0.39550781     ok        

# Example 4: Filter calculations with a reference list in a chosen group
candidates <- c("N", "O", "P", "Q", "T", "U", "V", "X", "Y")
rel_subset <- pedrel(tp, by = "Gen", reference = candidates)
#> Warning: pedrel(): 2 of 6 groups returned a non-ok status. Inspect the 'Status' and 'Message' columns. First issue: Gen = 3; Group has less than 2 individuals after applying 'reference'.
print(rel_subset)
#>      Gen NTotal NUsed  MeanRel  Status
#>    <int>  <int> <int>    <num>  <char>
#> 1:     1      9     2 0.000000      ok
#> 2:     2      5     2 0.500000      ok
#> 3:     3      7     1       NA skipped
#> 4:     4      3     2 0.156250      ok
#> 5:     5      2     2 0.046875      ok
#> 6:     6      2     0       NA skipped
#>                                                          Message
#>                                                           <char>
#> 1:                                                              
#> 2:                                                              
#> 3: Group has less than 2 individuals after applying 'reference'.
#> 4:                                                              
#> 5:                                                              
#> 6: Group has less than 2 individuals after applying 'reference'.
# }
```
