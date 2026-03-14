# 4. Pedigree Analysis and Population Genetics

This vignette summarizes a practical workflow for pedigree analysis with
`visPedigree`, with emphasis on the interpretation of the main
indicators used in breeding and conservation genetics.

The discussion is organized around five questions:

1.  How complete and deep is the pedigree?
2.  How long is the generation cycle?
3.  Has genetic diversity been eroded by unequal founder use,
    bottlenecks, or drift?
4.  How large is the effective population size under different
    definitions?
5.  Are relationship, inbreeding, subpopulation structure, and gene flow
    under control?

## 1. Setup and Data Preparation

Different package datasets are used for different analytical tasks.

- `deep_ped` is useful for pedigree depth and diversity summaries.
- `big_family_size_ped` contains a `Year` column and is suitable for
  generation intervals.
- `small_ped` is convenient for relationship examples.
- `inbred_ped` is convenient for inbreeding and classification examples.

``` r
library(visPedigree)
library(data.table)

data(deep_ped, package = "visPedigree")
data(big_family_size_ped, package = "visPedigree")
data(small_ped, package = "visPedigree")
data(inbred_ped, package = "visPedigree")

tp_deep <- tidyped(deep_ped)
tp_small <- tidyped(small_ped)
tp_inbred <- tidyped(inbred_ped)
```

## 2. Pedigree Overview with `pedstats()`

[`pedstats()`](https://luansheng.github.io/visPedigree/reference/pedstats.md)
provides a compact structural summary with three components:

- `$summary`: pedigree size and parental structure.
- `$ecg`: pedigree completeness and ancestral depth.
- `$gen_intervals`: generation intervals, only when a usable time column
  is available.

Here `deep_ped` has no explicit birth-date column. Accordingly,
[`pedstats()`](https://luansheng.github.io/visPedigree/reference/pedstats.md)
returns `$summary` and `$ecg`, whereas `$gen_intervals` remains `NULL`.

``` r
stats_deep <- pedstats(tp_deep)

stats_deep$summary
#>        N NSire  NDam NFounder MaxGen
#>    <int> <int> <int>    <int>  <int>
#> 1:  4399   483   554      138     13
tail(stats_deep$ecg)
#>         Ind      ECG FullGen MaxGen
#>      <char>    <num>   <num>  <num>
#> 1: K110997Q 7.616211       5     12
#> 2: K110997Z 6.722656       5     12
#> 3: K110998Q 4.606445       2     12
#> 4: K110998Z 6.417969       5     12
#> 5: K110999Q 7.345215       5     12
#> 6: K110999Z 6.875977       5     12
stats_deep$gen_intervals
#> NULL
```

Pedigree depth and pedigree time are related but distinct quantities.
The former is addressed by
[`pedecg()`](https://luansheng.github.io/visPedigree/reference/pedecg.md),
whereas the latter is addressed by
[`pedgenint()`](https://luansheng.github.io/visPedigree/reference/pedgenint.md).

## 3. Pedigree Completeness with `pedecg()`

Equivalent complete generations (ECG) summarize the amount of ancestral
information available for each individual:

$$ECG_{i} = \sum\limits_{j}\left( \frac{1}{2} \right)^{g_{ij}}$$

where $g_{ij}$ is the number of generations between individual $i$ and
ancestor $j$. ECG increases with both pedigree depth and pedigree
completeness.

In practice:

- `ECG` combines depth and completeness.
- `FullGen` counts how many fully known ancestral generations exist.
- `MaxGen` records the deepest known ancestral path.

`ECG` is especially useful because it provides the depth adjustment
required by realized effective population size estimators based on
inbreeding or coancestry.

``` r
ecg_deep <- pedecg(tp_deep)

ecg_deep[order(-ECG)][1:10]
#>          Ind      ECG FullGen MaxGen
#>       <char>    <num>   <num>  <num>
#>  1: K110034Q 9.078125       6     12
#>  2: K110052L 9.078125       6     12
#>  3: K110060H 9.078125       6     12
#>  4: K110069Z 9.078125       6     12
#>  5: K110097Q 9.078125       6     12
#>  6: K110118M 9.078125       6     12
#>  7: K110131M 9.078125       6     12
#>  8: K110138M 9.078125       6     12
#>  9: K110155M 9.078125       6     12
#> 10: K110165M 9.078125       6     12
```

## 4. Generation Intervals with `pedgenint()`

Generation interval is the age of a parent at the birth of its
offspring:

$$L = t_{offspring} - t_{parent}$$

[`pedgenint()`](https://luansheng.github.io/visPedigree/reference/pedgenint.md)
estimates this quantity for seven pathway summaries:

- `SS`, `SD`, `DS`, `DD`: sex-specific gametic pathways.
- `SO`, `DO`: sex-independent sire-offspring and dam-offspring pathways.
- `Average`: all parent-offspring pairs combined.

The function accepts `Date`/`POSIXct` columns, date strings, or numeric
years. When only an integer year is available, it is converted
internally to `YYYY-07-01` as a mid-year approximation.

``` r
tp_time <- tidyped(big_family_size_ped)

genint_year <- pedgenint(tp_time, timevar = "Year", unit = "year")
#> Numeric time column detected. Converting to Date (YYYY-07-01). For finer precision, convert to Date beforehand.
genint_year
#>    Pathway      N     Mean         SD
#>     <char>  <int>    <num>      <num>
#> 1: Average 280512 1.001093 0.03770707
#> 2:      DD    607 1.164398 0.37080261
#> 3:      DO 140256 1.001093 0.03770714
#> 4:      DS    507 1.196959 0.39776787
#> 5:      SD    607 1.164398 0.37080261
#> 6:      SO 140256 1.001093 0.03770714
#> 7:      SS    507 1.196959 0.39776787
```

The optional `cycle` parameter adds `GenEquiv`, which compares the
observed mean interval with a target breeding cycle:

$$GenEquiv_{i} = \frac{{\bar{L}}_{i}}{L_{cycle}}$$

Values larger than 1 indicate that the realized generation interval
exceeds the target cycle.

``` r
genint_cycle <- pedgenint(tp_time, timevar = "Year", unit = "year", cycle = 1.2)
#> Numeric time column detected. Converting to Date (YYYY-07-01). For finer precision, convert to Date beforehand.

genint_cycle[Pathway %in% c("SS", "SD", "DS", "DD", "Average")]
#>    Pathway      N     Mean         SD  GenEquiv
#>     <char>  <int>    <num>      <num>     <num>
#> 1: Average 280512 1.001093 0.03770707 0.8342445
#> 2:      DD    607 1.164398 0.37080261 0.9703321
#> 3:      DS    507 1.196959 0.39776787 0.9974660
#> 4:      SD    607 1.164398 0.37080261 0.9703321
#> 5:      SS    507 1.196959 0.39776787 0.9974660
```

## 5. Subpopulation Structure with `pedsubpop()`

Before interpreting diversity or relationship metrics, it is useful to
check whether the pedigree forms a single connected population or a
mixture of separate components.

[`pedsubpop()`](https://luansheng.github.io/visPedigree/reference/pedsubpop.md)
has two modes:

- `by = NULL`: summarize disconnected pedigree components.
- `by = "..."`: summarize the pedigree by an existing grouping variable.

``` r
ped_demo <- data.table(
  Ind = c("A", "B", "C", "D", "E", "F", "G", "H"),
  Sire = c(NA, NA, "A", NA, NA, "E", NA, NA),
  Dam  = c(NA, NA, "B", NA, NA, NA, NA, NA),
  Sex = c("male", "female", "male", "female", "male", "female", "male", "female"),
  Batch = c("L1", "L1", "L1", "L1", "L2", "L2", "L3", "L3")
)

tp_demo <- tidyped(ped_demo)

pedsubpop(tp_demo)
#>       Group     N N_Sire N_Dam N_Founder
#>      <char> <int>  <num> <num>     <int>
#> 1:      GP1     3      1     1         2
#> 2:      GP2     2      1     0         1
#> 3: Isolated     3      0     0         3
pedsubpop(tp_demo, by = "Batch")
#>     Group     N N_Sire N_Dam N_Founder
#>    <char> <int>  <int> <int>     <int>
#> 1:     L1     4      1     1         3
#> 2:     L3     2      0     0         2
#> 3:     L2     2      1     0         1
```

This is a compact summary tool. When the actual split pedigree objects
are needed for downstream analysis, use
[`splitped()`](https://luansheng.github.io/visPedigree/reference/splitped.md).

## 6. Diversity Indicators with `pediv()`

[`pediv()`](https://luansheng.github.io/visPedigree/reference/pediv.md)
is the integrated diversity summary. It combines:

- founder and ancestor contributions (`fe`, `fa`),
- founder genome equivalents (`fg`),
- three effective population size estimates (`Ne`).

All of these quantities depend on the definition of the **reference
population**. In the present example, the reference population is
defined as the most recent two generations in `deep_ped`.

``` r
ref_pop <- tp_deep[Gen >= max(Gen) - 1, Ind]
length(ref_pop)
#> [1] 3471
```

``` r
div_res <- pediv(tp_deep, reference = ref_pop, top = 10, seed = 42L)
#> Calculating founder and ancestor contributions...
#> Calculating founder contributions...
#> Calculating ancestor contributions (Boichard's iterative algorithm)...
#> Calculating Ne (coancestry) and fg...
#> Calculating Ne (inbreeding)...
#> Calculating Ne (demographic)...

div_res$summary
#>     NRef NFounder       fe NAncestor       fa     fafe       fg  MeanCoan
#>    <int>    <int>    <num>     <int>    <num>    <num>    <num>     <num>
#> 1:  3471      157 64.73344        94 44.12033 0.681569 19.17965 0.0260693
#>    NSampledCoan NeCoancestry NeInbreeding NeDemographic
#>           <int>        <num>        <num>         <num>
#> 1:         1000     124.5124     98.00425      374.2021
div_res$ancestors
#>          Ind    Contrib CumContrib  Rank
#>       <char>      <num>      <num> <int>
#>  1: K500I804 0.05248848 0.05248848     1
#>  2: K60GXQ91 0.04525893 0.09774741     2
#>  3: K60NXQ91 0.04525893 0.14300634     3
#>  4: K600532I 0.04445765 0.18746399     4
#>  5: K700S069 0.03927182 0.22673581     5
#>  6: K900D251 0.03622875 0.26296456     6
#>  7: K700S011 0.02906223 0.29202679     7
#>  8: K700U416 0.02890017 0.32092697     8
#>  9: K700U650 0.02890017 0.34982714     9
#> 10: K500I869 0.02831497 0.37814211    10
```

### 6.1 `fe`, `fa`, and `fg`: what do they measure?

These three indicators describe diversity loss from complementary
angles.

#### Effective number of founders (`fe`)

$$f_{e} = \frac{1}{\sum\limits_{i = 1}^{f}p_{i}^{2}}$$

where $p_{i}$ is the expected contribution of founder $i$ to the
reference population. `fe` answers the question: if all founders had
contributed equally, how many founders would generate the same
diversity?

Use `fe` when the main concern is **unequal founder representation**.

#### Effective number of ancestors (`fa`)

$$f_{a} = \frac{1}{\sum\limits_{j = 1}^{a}q_{j}^{2}}$$

where $q_{j}$ is the **marginal** contribution of ancestor $j$ after
removing the contributions already explained by more influential
ancestors. `fa` is usually smaller than `fe` when bottlenecks occurred.

Use `fa` when the main concern is **genetic bottlenecks caused by a
limited set of influential ancestors**.

#### Founder genome equivalents (`fg`)

In its simplest interpretation,

$$f_{g} = \frac{1}{2\bar{C}}$$

where $\bar{C}$ is the mean coancestry of the reference population. In
`visPedigree`, `fg` is estimated from a diagonal-corrected mean
coancestry:

$$\widehat{\bar{C}} = \frac{N - 1}{N} \cdot \frac{{\bar{a}}_{off}}{2} + \frac{1 + {\bar{F}}_{s}}{2N}$$

$$f_{g} = \frac{1}{2\widehat{\bar{C}}}$$

where $N$ is the size of the full reference cohort, ${\bar{a}}_{off}$ is
the mean off-diagonal additive relationship among sampled individuals,
and ${\bar{F}}_{s}$ is their mean inbreeding coefficient.

Use `fg` when the main concern is **the amount of founder genome still
surviving after unequal use, bottlenecks, and drift**. In practice, `fg`
is often the most conservative diversity indicator.

## 7. Effective Population Size with `pedne()` and `pediv()`

[`pediv()`](https://luansheng.github.io/visPedigree/reference/pediv.md)
reports three complementary effective population size definitions. Each
addresses a distinct biological question.

### 7.1 Demographic `Ne`

$$N_{e} = \frac{4N_{m}N_{f}}{N_{m} + N_{f}}$$

where $N_{m}$ and $N_{f}$ are the numbers of contributing males and
females. This is the easiest estimate to understand, but it ignores
realized pedigree structure, inbreeding, and drift. It is therefore
often optimistic.

### 7.2 Inbreeding `Ne`

$$\Delta F_{i} = 1 - \left( 1 - F_{i} \right)^{1/{(ECG_{i} - 1)}}$$

$$N_{e} = \frac{1}{2\overline{\Delta F}}$$

This estimator uses the realized rate of inbreeding. `ECG` standardizes
individuals with different pedigree depths. Use this estimate when the
primary concern is **the rate of inbreeding accumulation**.

### 7.3 Coancestry `Ne`

$$\Delta c_{ij} = 1 - \left( 1 - c_{ij} \right)^{1/{(\frac{ECG_{i} + ECG_{j}}{2})}}$$

$$N_{e} = \frac{1}{2\overline{\Delta c}}$$

This estimator is based on the rate of coancestry among members of the
reference population. Because it captures the accumulation of
relatedness before it is fully expressed as realized inbreeding, it is
often the strictest and most sensitive warning signal in managed
breeding populations.

## 8. Average Relationship Trends with `pedrel()`

[`pedrel()`](https://luansheng.github.io/visPedigree/reference/pedrel.md)
summarizes the mean off-diagonal additive relationship within groups:

$$MeanRel = \frac{1}{n(n - 1)}\sum\limits_{i \neq j}a_{ij}$$

where $a_{ij} = 2f_{ij}$. A higher `MeanRel` means that individuals
within the group are, on average, more related by descent.

This summary is useful for tracking relatedness by generation or year.

``` r
tp_small$BirthYear <- 2010 + tp_small$Gen

rel_by_gen <- pedrel(tp_small, by = "Gen")
rel_by_gen
#>      Gen NTotal NUsed   MeanRel
#>    <int>  <int> <int>     <num>
#> 1:     1      9     9 0.0000000
#> 2:     2      5     5 0.2000000
#> 3:     3      7     7 0.1547619
#> 4:     4      3     3 0.1145833
#> 5:     5      2     2 0.0468750
#> 6:     6      2     2 0.5507812
```

The `reference` argument lets you focus on a subset of interest inside
each group, such as candidate breeders.

``` r
ref_ids <- c("Z1", "Z2", "X", "Y")

pedrel(tp_small, by = "Gen", reference = ref_ids)
#> Warning in FUN(X[[i]], ...): Group '1' has less than 2 individuals after
#> applying 'reference', returning NA_real_.
#> Warning in FUN(X[[i]], ...): Group '2' has less than 2 individuals after
#> applying 'reference', returning NA_real_.
#> Warning in FUN(X[[i]], ...): Group '3' has less than 2 individuals after
#> applying 'reference', returning NA_real_.
#> Warning in FUN(X[[i]], ...): Group '4' has less than 2 individuals after
#> applying 'reference', returning NA_real_.
#>      Gen NTotal NUsed   MeanRel
#>    <int>  <int> <int>     <num>
#> 1:     1      9     0        NA
#> 2:     2      5     0        NA
#> 3:     3      7     0        NA
#> 4:     4      3     0        NA
#> 5:     5      2     2 0.0468750
#> 6:     6      2     2 0.5507812
```

`MeanRel` and coancestry-based `Ne` are conceptually linked, but they
are not identical summaries. `MeanRel` is an average additive
relationship within a group, whereas coancestry-based `Ne` is derived
from the *rate of increase* in coancestry across pedigree depth.

## 9. Inbreeding Trends with `inbreed()` and `pedfclass()`

Wright’s inbreeding coefficient $F$ is the probability that the two
alleles of an individual are identical by descent. A practical starting
point is to inspect the mean trend by generation.

``` r
tp_inbred_f <- inbreed(tp_inbred)

f_trend <- tp_inbred_f[, .(MeanF = mean(f, na.rm = TRUE)), by = Gen]
f_trend
#> Tidy Pedigree Object
#>      Gen  MeanF
#>    <int>  <num>
#> 1:     1 0.0000
#> 2:     2 0.0000
#> 3:     3 0.2500
#> 4:     4 0.2500
#> 5:     5 0.4375
```

For classification by inbreeding severity,
[`pedfclass()`](https://luansheng.github.io/visPedigree/reference/pedfclass.md)
can be applied directly to the pedigree object. If the inbreeding
coefficient is not yet present, the function computes it internally.

``` r
pedfclass(tp_inbred)
#> Calculating inbreeding coefficients...
#> Key: <FClass>
#>                 FClass Count Percentage
#>                  <ord> <int>      <num>
#> 1:               F = 0     4   57.14286
#> 2:     0 < F <= 0.0625     0    0.00000
#> 3: 0.0625 < F <= 0.125     0    0.00000
#> 4:   0.125 < F <= 0.25     2   28.57143
#> 5:            F > 0.25     1   14.28571
```

Custom reporting thresholds can be supplied through `breaks`.

``` r
pedfclass(tp_inbred, breaks = c(0.03125, 0.0625, 0.125, 0.25))
#> Calculating inbreeding coefficients...
#> Key: <FClass>
#>                   FClass Count Percentage
#>                    <ord> <int>      <num>
#> 1:                 F = 0     4   57.14286
#> 2:      0 < F <= 0.03125     0    0.00000
#> 3: 0.03125 < F <= 0.0625     0    0.00000
#> 4:   0.0625 < F <= 0.125     0    0.00000
#> 5:     0.125 < F <= 0.25     2   28.57143
#> 6:              F > 0.25     1   14.28571
```

The default thresholds correspond approximately to familiar pedigree
scenarios:

- $F = 0.0625$: half-sib mating.
- $F = 0.125$: avuncular or grandparent-grandchild mating.
- $F = 0.25$: full-sib or parent-offspring mating.

## 10. Gene Flow and Partial Inbreeding

### 10.1 `pedancestry()`: founder-line proportions

[`pedancestry()`](https://luansheng.github.io/visPedigree/reference/pedancestry.md)
tracks how founder groups contribute to later descendants. This is
useful when founders are labeled by line, strain, or geographic origin.

``` r
ped_line <- data.table(
  Ind = c("A", "B", "C", "D", "E", "F", "G"),
  Sire = c(NA, NA, NA, NA, "A", "C", "E"),
  Dam  = c(NA, NA, NA, NA, "B", "D", "F"),
  Sex = c("male", "female", "male", "female", "male", "female", "male"),
  Line = c("Line1", "Line1", "Line2", "Line2", NA, NA, NA)
)

tp_line <- tidyped(ped_line)

anc <- pedancestry(tp_line, foundervar = "Line")
anc
#>       Ind Line1 Line2
#>    <char> <num> <num>
#> 1:      A   1.0   0.0
#> 2:      B   1.0   0.0
#> 3:      C   0.0   1.0
#> 4:      D   0.0   1.0
#> 5:      E   1.0   0.0
#> 6:      F   0.0   1.0
#> 7:      G   0.5   0.5
```

### 10.2 `pedpartial()`: which ancestors explain inbreeding?

[`pedpartial()`](https://luansheng.github.io/visPedigree/reference/pedpartial.md)
decomposes total inbreeding into contributions from selected ancestors.
It is useful for identifying which ancestors are most responsible for
the observed inbreeding.

``` r
partial_f <- pedpartial(tp_inbred, ancestors = c("A", "B"))
#> Calculating partial inbreeding for 2 ancestors...
partial_f
#>       Ind       A       B
#>    <char>   <num>   <num>
#> 1:      A 0.00000 0.00000
#> 2:      B 0.00000 0.00000
#> 3:      C 0.00000 0.00000
#> 4:      D 0.00000 0.00000
#> 5:      E 0.12500 0.12500
#> 6:      F 0.00000 0.25000
#> 7:      G 0.15625 0.28125
```

## 11. Practical Interpretation

One useful interpretation order is:

1.  Check pedigree depth and completeness with
    [`pedecg()`](https://luansheng.github.io/visPedigree/reference/pedecg.md).
2.  Check generation timing with
    [`pedgenint()`](https://luansheng.github.io/visPedigree/reference/pedgenint.md).
3.  Quantify diversity loss with `fe`, `fa`, and `fg`.
4.  Compare demographic, inbreeding, and coancestry-based `Ne`.
5.  Monitor `MeanRel` and `MeanF` over time.
6.  Use
    [`pedsubpop()`](https://luansheng.github.io/visPedigree/reference/pedsubpop.md),
    [`pedancestry()`](https://luansheng.github.io/visPedigree/reference/pedancestry.md),
    and
    [`pedpartial()`](https://luansheng.github.io/visPedigree/reference/pedpartial.md)
    to diagnose structure, introgression, and bottlenecks.

## References

- Boichard, D., Maignel, L., & Verrier, É. (1997). The value of using
  probabilities of gene origin to measure genetic variability in a
  population. *Genetics Selection Evolution*, 29(1), 5-23.
- Caballero, A., & Toro, M. A. (2000). Interrelations between effective
  population size and other pedigree tools for the management of
  conserved populations. *Genetical Research*, 75(3), 331-343.
- Cervantes, I., Goyache, F., Molina, A., Valera, M., & Gutiérrez, J. P.
  (2011). Estimation of effective population size from the rate of
  coancestry in pedigreed populations. *Journal of Animal Breeding and
  Genetics*, 128(1), 56-63.
- Gutiérrez, J. P., Cervantes, I., Molina, A., Valera, M., & Goyache, F.
  (2008). Individual increase in inbreeding allows estimating effective
  sizes from pedigrees. *Genetics Selection Evolution*, 40(4), 359-370.
- Gutiérrez, J. P., Cervantes, I., & Goyache, F. (2009). Improving the
  estimation of realized effective population sizes in farm animals.
  *Journal of Animal Breeding and Genetics*, 126(4), 327-332.
- Lacy, R. C. (1989). Analysis of founder representation in pedigrees:
  founder equivalents and founder genome equivalents. *Zoo Biology*,
  8(2), 111-123.
- Wright, S. (1922). Coefficients of inbreeding and relationship. *The
  American Naturalist*, 56(645), 330-338.
- Wright, S. (1931). Evolution in Mendelian populations. *Genetics*,
  16(2), 97-159.
