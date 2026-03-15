# Calculate Effective Population Size

Calculates the effective population size (Ne) based on the rate of
coancestry, the rate of inbreeding, or demographic parent numbers.

## Usage

``` r
pedne(
  ped,
  method = c("coancestry", "inbreeding", "demographic"),
  by = NULL,
  reference = NULL,
  nsamples = 1000,
  ncores = 1,
  seed = NULL
)
```

## Arguments

- ped:

  A `tidyped` object.

- method:

  Character. The method to compute Ne. One of `"coancestry"` (default),
  `"inbreeding"`, or `"demographic"`.

- by:

  Character. The name of the column used to group cohorts (e.g., "Year",
  "BirthYear"). If NULL, calculates overall Ne for all individuals.

- reference:

  Character vector. Optional subset of individual IDs defining the
  reference cohort. If NULL, uses all individuals in the pedigree.

- nsamples:

  Integer. Number of individuals to randomly sample per cohort when
  using the `"coancestry"` method. Very large cohorts will be sampled
  down to this size to save memory and time (default: 1000).

- ncores:

  Integer. Number of cores for parallel processing. Currently only
  effective for `method = "coancestry"` (default: 1).

- seed:

  Integer or NULL. Random seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before sampling in
  the coancestry method, ensuring reproducible \\N_e\\ estimates.
  Default is `NULL` (no fixed seed).

## Value

A `data.table` with columns:

- `Cohort`: Cohort or grouping variable value.

- `N`: Number of individuals in the cohort.

- `Ne`: Calculated effective population size.

- `...`: Additional columns depending on the selected method (e.g.,
  `NSampled`, `DeltaC`, `MeanF`, `DeltaF`, `Nm`, `Nf`).

## Details

The effective population size can be calculated using one of three
methods:

- **"coancestry"** (Default): Based on the rate of coancestry
  (\\f\_{ij}\\) between pairs of individuals. In this context,
  \\f\_{ij}\\ is the probability that two alleles randomly sampled from
  individuals \$i\$ and \$j\$ are identical by descent, which is
  equivalent to half the additive genetic relationship (\\f\_{ij} =
  a\_{ij} / 2\\). This method is generally more robust as it accounts
  for full genetic drift and bottlenecks (Cervantes et al., 2011).
  \$\$\Delta c\_{ij} = 1 - (1 - c\_{ij})^{1/(\frac{ECG_i +
  ECG_j}{2})}\$\$ \$\$N_e = \frac{1}{2 \overline{\Delta c}}\$\$ To
  handle large populations, this method samples `nsamples` individuals
  per cohort and computes the mean rate of coancestry among them.

- **"inbreeding"**: Based on the individual rate of inbreeding (\$F_i\$)
  (Gutiérrez et al., 2008, 2009). \$\$\Delta F_i = 1 - (1 -
  F_i)^{1/(ECG_i - 1)}\$\$ \$\$N_e = \frac{1}{2 \overline{\Delta F}}\$\$

- **"demographic"**: Based on the demographic census of breeding males
  and females (Wright, 1931). \$\$N_e = \frac{4 N_m N_f}{N_m + N_f}\$\$
  Where \\N_m\\ and \\N_f\\ are the number of unique male and female
  parents contributing to the cohort.

## References

Cervantes, I., Goyache, F., Molina, A., Valera, M., & Gutiérrez, J. P.
(2011). Estimation of effective population size from the rate of
coancestry in pedigreed populations. *Journal of Animal Breeding and
Genetics*, 128(1), 56-63.

Gutiérrez, J. P., Cervantes, I., Molina, A., Valera, M., & Goyache, F.
(2008). Individual increase in inbreeding allows estimating effective
sizes from pedigrees. *Genetics Selection Evolution*, 40(4), 359-370.

Gutiérrez, J. P., Cervantes, I., & Goyache, F. (2009). Improving the
estimation of realized effective population sizes in farm animals.
*Journal of Animal Breeding and Genetics*, 126(4), 327-332.

Wright, S. (1931). Evolution in Mendelian populations. *Genetics*,
16(2), 97-159.

## Examples

``` r
# \donttest{
# Coancestry-based Ne (default) using a simple pedigree grouped by year
tp_simple <- tidyped(simple_ped)
tp_simple$BirthYear <- 2000 + tp_simple$Gen
ne_coan <- suppressMessages(pedne(tp_simple, by = "BirthYear", seed = 42L))
ne_coan
#>    Cohort     N NSampled       DeltaC       Ne
#>     <num> <int>    <int>        <num>    <num>
#> 1:   2001    28       28           NA       NA
#> 2:   2002    16       16 0.0031250000 160.0000
#> 3:   2003     8        8 0.0017708540 282.3496
#> 4:   2004     4        4 0.0009737864 513.4596
#> 5:   2005     2        2 0.0011276950 443.3823
#> 6:   2006     1        1           NA       NA

# Inbreeding-based Ne using an inbred pedigree
tp_inbred <- tidyped(inbred_ped)
ne_inb <- suppressMessages(pedne(tp_inbred, method = "inbreeding", by = "Gen"))
ne_inb
#> Tidy Pedigree Object
#>    Cohort     N  MeanF DeltaF MeanECG    Ne
#>     <int> <int>  <num>  <num>   <num> <num>
#> 1:      3     1 0.2500   0.25       2     2
#> 2:      4     1 0.2500   0.25       2     2
#> 3:      5     1 0.4375   0.25       3     2

# Demographic Ne from the number of contributing sires and dams
ne_demo <- suppressMessages(pedne(tp_simple, method = "demographic", by = "BirthYear"))
ne_demo
#> Tidy Pedigree Object
#>    Cohort     N    Nm    Nf        Ne
#>     <num> <int> <int> <int>     <num>
#> 1:   2001    28     0     0        NA
#> 2:   2002    16    14    14 28.000000
#> 3:   2003     8     7     8 14.933333
#> 4:   2004     4     4     3  6.857143
#> 5:   2005     2     2     2  4.000000
#> 6:   2006     1     1     1  2.000000
# }
```
