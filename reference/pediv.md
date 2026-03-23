# Calculate Genetic Diversity Indicators

Combines founder/ancestor contributions (\$f_e\$, \$f_a\$) and effective
population size estimates (Ne) from three methods into a single summary
object.

## Usage

``` r
pediv(
  ped,
  reference = NULL,
  top = 20,
  nsamples = 1000,
  ncores = 1,
  seed = NULL
)
```

## Arguments

- ped:

  A `tidyped` object.

- reference:

  Character vector. Optional subset of individual IDs defining the
  reference population. If NULL, uses all individuals in the most recent
  generation.

- top:

  Integer. Number of top contributors to return in founder/ancestor
  tables. Default is 20.

- nsamples:

  Integer. Number of individuals sampled per cohort for the coancestry
  Ne method and for \\f_g\\ estimation. Very large cohorts are sampled
  down to this size to control memory usage (default: 1000).

- ncores:

  Integer. Number of cores for parallel processing in the coancestry
  method. Default is 1.

- seed:

  Integer or NULL. Random seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before sampling in
  the coancestry method, ensuring reproducible \\f_g\\ and \\N_e\\
  estimates. Default is `NULL` (no fixed seed).

## Value

A list with class `pediv` containing:

- `summary`: A single-row `data.table` with columns `NRef`, `NFounder`,
  `feH`, `fe`, `NAncestor`, `faH`, `fa`, `fafe`, `fg`, `MeanCoan`,
  `GeneDiv`, `NSampledCoan`, `NeCoancestry`, `NeInbreeding`,
  `NeDemographic`. Here `feH` and `faH` are the Shannon-entropy
  (\\q=1\\) effective numbers of founders and ancestors, respectively.
  `GeneDiv` is the pedigree-based retained genetic diversity, computed
  as \\1 - \bar{C}\\, where \\\bar{C}\\ is the diagonal-corrected
  population mean coancestry (`MeanCoan`).

- `founders`: A `data.table` of top founder contributions.

- `ancestors`: A `data.table` of top ancestor contributions.

## Details

Internally calls
[`pedcontrib`](https://luansheng.github.io/visPedigree/reference/pedcontrib.md)
for \\f_e\\ and \\f_a\\. The coancestry method is called via the
internal `calc_ne_coancestry()` function directly so that \\f_g\\ and
the Ne estimate can be obtained from the same traced pedigree without
duplication. The inbreeding and demographic Ne methods are obtained via
[`pedne`](https://luansheng.github.io/visPedigree/reference/pedne.md).
All calculations use the same `reference` population. If any method
fails (e.g., insufficient pedigree depth), its value is `NA` rather than
stopping execution.

\\f_g\\ (founder genome equivalents, Caballero & Toro 2000) is estimated
from the diagonal-corrected mean coancestry of the reference population:
\$\$\hat{\bar{C}} = \frac{N-1}{N} \cdot \frac{\bar{a}\_{off}}{2} +
\frac{1 + \bar{F}\_s}{2N}\$\$ \$\$f_g = \frac{1}{2 \hat{\bar{C}}}\$\$
where \\N\\ is the full reference cohort size, \\\bar{a}\_{off}\\ is the
off-diagonal mean relationship among sampled individuals, and
\\\bar{F}\_s\\ is their mean inbreeding coefficient.

## See also

[`pedcontrib`](https://luansheng.github.io/visPedigree/reference/pedcontrib.md),
[`pedne`](https://luansheng.github.io/visPedigree/reference/pedne.md),
[`pedstats`](https://luansheng.github.io/visPedigree/reference/pedstats.md)

## Examples

``` r
# \donttest{
tp <- tidyped(small_ped)
div <- pediv(tp, reference = c("Z1", "Z2", "X", "Y"), seed = 42L)
#> Calculating founder and ancestor contributions...
#> Calculating founder contributions...
#> Calculating ancestor contributions (Boichard's iterative algorithm)...
#> Calculating Ne (coancestry) and fg...
#> Calculating Ne (inbreeding)...
#> Calculating Ne (demographic)...
print(div)
#> Genetic Diversity Summary
#> =========================
#> Reference population size : 4
#> 
#> -- Founder / Ancestor Contributions --
#> Founders  (total) : 9    fe(H) = 7.672    fe = 6.585
#> Ancestors (total) : 3    fa(H) = 2.828    fa = 2.667
#> fa/fe ratio       : 0.4049
#> Founder genomes   : fg = 2.075  (MeanCoan = 0.240967, NSampled = 4)
#> Gene diversity    : GeneDiv = 0.7590  (= 1 - MeanCoan)
#> Hierarchy: fg <= fa <= fe <= NFounder  =  2.07 <= 2.667 <= 6.585 <= 9
#> 
#> -- Effective Population Size (Ne) --
#>   Coancestry  :     9.0
#>   Inbreeding  :    29.4
#>   Demographic :     6.0
#> 
#> Top Founder Contributions (top 5):
#>       Ind  Contrib CumContrib  Rank
#>    <char>    <num>      <num> <int>
#> 1:      N 0.281250   0.281250     1
#> 2:     J2 0.125000   0.406250     2
#> 3:      R 0.125000   0.531250     3
#> 4:      A 0.109375   0.640625     4
#> 5:      B 0.109375   0.750000     5
#> 
#> Top Ancestor Contributions (top 5):
#>       Ind Contrib CumContrib  Rank
#>    <char>   <num>      <num> <int>
#> 1:      X    0.50       0.50     1
#> 2:      N    0.25       0.75     2
#> 3:      Y    0.25       1.00     3

# Access Shannon effective numbers from summary
div$summary$feH   # Shannon effective founders (q=1)
#> [1] 7.671504
div$summary$faH   # Shannon effective ancestors (q=1)
#> [1] 2.828427

# Founder diversity profile: NFounder >= feH >= fe
with(div$summary, c(NFounder = NFounder, feH = feH, fe = fe))
#> NFounder      feH       fe 
#> 9.000000 7.671504 6.585209 
# }
```
