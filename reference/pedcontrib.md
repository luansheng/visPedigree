# Calculate Founder and Ancestor Contributions

Calculates genetic contributions from founders and influential
ancestors. Implements the gene dropping algorithm for founder
contributions and Boichard's algorithm for ancestor contributions to
estimate the effective number of founders (\$f_e\$) and ancestors
(\$f_a\$).

## Usage

``` r
pedcontrib(
  ped,
  reference = NULL,
  mode = c("both", "founder", "ancestor"),
  top = 20
)
```

## Arguments

- ped:

  A `tidyped` object.

- reference:

  Character vector. Optional subset of individual IDs defining the
  reference population. If NULL, uses all individuals in the most recent
  generation.

- mode:

  Character. Type of contribution to calculate:

  - `"founder"`: Founder contributions (\$f_e\$).

  - `"ancestor"`: Ancestor contributions (\$f_a\$).

  - `"both"`: Both founder and ancestor contributions.

- top:

  Integer. Number of top contributors to return. Default is 20.

## Value

A list with class `pedcontrib` containing:

- `founders`: A `data.table` of founder contributions (if mode includes
  "founder", or "both").

- `ancestors`: A `data.table` of ancestor contributions (if mode
  includes "ancestor", or "both").

- `summary`: A `list` of summary statistics including the effective
  number of founders (\$f_e\$) and ancestors (\$f_a\$).

Each contribution table contains:

- `Ind`: Individual ID.

- `Contrib`: Contribution to the reference population (0-1).

- `CumContrib`: Cumulative contribution.

- `Rank`: Rank by contribution.

## Details

\*\*Founder Contributions (\$f_e\$):\*\* Calculated by probabilistic
gene flow from founders to the reference cohort. When individual
ancestors with one unknown parent exist, "phantom" parents are
temporarily injected correctly conserving the probability mass.

\*\*Ancestor Contributions (\$f_a\$):\*\* Calculated using Boichard's
iterative algorithm (1997), accounting for:

- Marginal genetic contribution of each ancestor

- Long-term contributions through multiple pathways

The parameter \$f_a\$ acts as a stringent metric since it identifies the
bottlenecks of genetic variation in pedigrees.

## References

Boichard, D., Maignel, L., & Verrier, É. (1997). The value of using
probabilities of gene origin to measure genetic variability in a
population. Genetics Selection Evolution, 29(1), 5-23.

## Examples

``` r
# \donttest{
library(data.table)
# Load a sample pedigree
tp <- tidyped(small_ped)

# Calculate both founder and ancestor contributions for reference population
ref_ids <- c("Z1", "Z2", "X", "Y")
contrib <- pedcontrib(tp, reference = ref_ids, mode = "both")
#> Calculating founder contributions...
#> Calculating ancestor contributions (Boichard's iterative algorithm)...

# Print results including f_e and f_a
print(contrib)
#> Founder and Ancestor Contributions
#> ===================================
#> Reference population size: 4
#> 
#> Founders: 9 (reported top 9)
#> Effective number of founders (f_e): 6.59
#> 
#> Top 10 Founder Contributions:
#>       Ind  Contrib CumContrib  Rank
#>    <char>    <num>      <num> <int>
#> 1:      N 0.281250   0.281250     1
#> 2:     J2 0.125000   0.406250     2
#> 3:      R 0.125000   0.531250     3
#> 4:      A 0.109375   0.640625     4
#> 5:      B 0.109375   0.750000     5
#> 6:      F 0.093750   0.843750     6
#> 7:      I 0.062500   0.906250     7
#> 8:     J1 0.062500   0.968750     8
#> 9:      O 0.031250   1.000000     9
#> 
#> Ancestors: 3 (reported top 3)
#> Effective number of ancestors (f_a): 2.67
#> 
#> Top 10 Ancestor Contributions:
#>       Ind Contrib CumContrib  Rank
#>    <char>   <num>      <num> <int>
#> 1:      X    0.50       0.50     1
#> 2:      N    0.25       0.75     2
#> 3:      Y    0.25       1.00     3
# }
```
