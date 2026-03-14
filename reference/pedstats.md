# Pedigree Statistics

Calculates comprehensive statistics for a pedigree, including population
structure, generation intervals, and ancestral depth.

## Usage

``` r
pedstats(
  ped,
  timevar = NULL,
  unit = "year",
  cycle = NULL,
  ecg = TRUE,
  genint = TRUE,
  ...
)
```

## Arguments

- ped:

  A `tidyped` object.

- timevar:

  Optional character. Name of the column containing the **birth date**
  (or hatch date) of each individual. Accepted column formats:

  - `Date` or `POSIXct` (recommended).

  - A date string parseable by `as.POSIXct` (e.g., `"2020-06-15"`). Use
    `format` via `...` for non-ISO strings.

  - A numeric year (e.g., `2020`). Automatically converted to `Date`
    (`"YYYY-07-01"`) with a message.

  If `NULL`, attempts auto-detection from common column names
  (`"BirthYear"`, `"Year"`, `"BirthDate"`, etc.).

- unit:

  Character. Time unit for reporting generation intervals: `"year"`
  (default), `"month"`, `"day"`, or `"hour"`.

- cycle:

  Numeric. Optional target generation cycle length in `unit`s. When
  provided, `gen_intervals` will include a `GenEquiv` column (observed
  Mean / cycle). See
  [`pedgenint`](https://luansheng.github.io/visPedigree/reference/pedgenint.md)
  for details.

- ecg:

  Logical. Whether to compute equivalent complete generations for each
  individual via
  [`pedecg`](https://luansheng.github.io/visPedigree/reference/pedecg.md).
  Default `TRUE`.

- genint:

  Logical. Whether to compute generation intervals via
  [`pedgenint`](https://luansheng.github.io/visPedigree/reference/pedgenint.md).
  Requires a detectable `timevar` column. Default `TRUE`.

- ...:

  Additional arguments passed to
  [`pedgenint`](https://luansheng.github.io/visPedigree/reference/pedgenint.md),
  e.g., `format` for custom date parsing or `by` for grouping.

## Value

An object of class `pedstats`, which is a list containing:

- `summary`: A `data.table` with one row summarising the whole pedigree.
  Columns:

  - `N` — total number of individuals.

  - `NSire` — number of unique sires.

  - `NDam` — number of unique dams.

  - `NFounder` — number of founder individuals (both parents unknown).

  - `MaxGen` — maximum generation number.

- `ecg`: A `data.table` with one row per individual (`NULL` if
  `ecg = FALSE`). Columns:

  - `Ind` — individual identifier.

  - `ECG` — equivalent complete generations.

  - `FullGen` — number of fully known generations.

  - `MaxGen` — maximum traceable generation depth.

- `gen_intervals`: A `data.table` of generation intervals (`NULL` if no
  `timevar` is detected or `genint = FALSE`). Columns:

  - `Pathway` — gametic pathway label. Seven values: `"SS"` (sire to
    son), `"SD"` (sire to daughter), `"DS"` (dam to son), `"DD"` (dam to
    daughter) — require offspring sex; `"SO"` (sire to offspring) and
    `"DO"` (dam to offspring) — sex-independent; and `"Average"` — all
    parent-offspring pairs combined.

  - `N` — number of parent-offspring pairs.

  - `Mean` — mean generation interval.

  - `SD` — standard deviation of the interval.

  - `GenEquiv` — `Mean / cycle` (only present when `cycle` is supplied).

## Examples

``` r
# \donttest{
# ---- Without time variable ----
tp <- tidyped(simple_ped)
ps <- pedstats(tp)
ps$summary
#>        N NSire  NDam NFounder MaxGen
#>    <int> <int> <int>    <int>  <int>
#> 1:    59    28    28       28      6
ps$ecg
#>        Ind     ECG FullGen MaxGen
#>     <char>   <num>   <num>  <num>
#>  1: J0C032 0.00000       0      0
#>  2: J0C185 0.00000       0      0
#>  3: J0C231 0.00000       0      0
#>  4: J0C317 0.00000       0      0
#>  5: J0C355 0.00000       0      0
#>  6: J0C450 0.00000       0      0
#>  7: J0C561 0.00000       0      0
#>  8: J0C583 0.00000       0      0
#>  9: J0C591 0.00000       0      0
#> 10: J0C612 0.00000       0      0
#> 11: J0Z060 0.00000       0      0
#> 12: J0Z167 0.00000       0      0
#> 13: J0Z256 0.00000       0      0
#> 14: J0Z333 0.00000       0      0
#> 15: J0Z380 0.00000       0      0
#> 16: J0Z444 0.00000       0      0
#> 17: J0Z475 0.00000       0      0
#> 18: J0Z482 0.00000       0      0
#> 19: J0Z511 0.00000       0      0
#> 20: J0Z563 0.00000       0      0
#> 21: J0Z624 0.00000       0      0
#> 22: J0Z664 0.00000       0      0
#> 23: J0Z808 0.00000       0      0
#> 24: J0Z839 0.00000       0      0
#> 25: J0Z843 0.00000       0      0
#> 26: J0Z848 0.00000       0      0
#> 27: J0Z938 0.00000       0      0
#> 28: J0Z990 0.00000       0      0
#> 29: J1C802 1.00000       1      1
#> 30: J1C929 1.00000       1      1
#> 31: J1E539 1.00000       1      1
#> 32: J1E852 1.00000       1      1
#> 33: J1F266 1.00000       1      1
#> 34: J1H419 1.00000       1      1
#> 35: J1H604 1.00000       1      1
#> 36: J1I438 1.00000       1      1
#> 37: J1I975 1.00000       1      1
#> 38: J1J134 1.00000       1      1
#> 39: J1J576 1.00000       1      1
#> 40: J1J858 0.50000       0      1
#> 41: J1K462 1.00000       1      1
#> 42: J1X971 1.00000       1      1
#> 43: J1Y339 1.00000       1      1
#> 44: J1Z417 1.00000       1      1
#> 45: J2C161 2.00000       2      2
#> 46: J2C808 2.00000       2      2
#> 47: J2F588 1.00000       0      2
#> 48: J2G465 2.00000       2      2
#> 49: J2X544 1.75000       1      2
#> 50: J2Y434 2.00000       2      2
#> 51: J2Z411 2.00000       2      2
#> 52: J2Z903 2.00000       2      2
#> 53: J3L886 2.50000       1      3
#> 54: J3X697 1.50000       0      3
#> 55: J3Y620 3.00000       3      3
#> 56: J3Y771 2.87500       2      3
#> 57: J4E185 3.00000       1      4
#> 58: J4Y326 3.93750       3      4
#> 59: J5X804 4.46875       2      5
#>        Ind     ECG FullGen MaxGen
#>     <char>   <num>   <num>  <num>

# ---- With annual Year column (big_family_size_ped) ----
tp2 <- tidyped(big_family_size_ped)
ps2 <- pedstats(tp2, timevar = "Year")
#> Numeric time column detected. Converting to Date (YYYY-07-01). For finer precision, convert to Date beforehand.
ps2$summary
#>         N NSire  NDam NFounder MaxGen
#>     <int> <int> <int>    <int>  <int>
#> 1: 178431   707   808      672      9
ps2$gen_intervals
#>    Pathway      N     Mean         SD
#>     <char>  <int>    <num>      <num>
#> 1: Average 280512 1.001093 0.03770707
#> 2:      DD    607 1.164398 0.37080261
#> 3:      DO 140256 1.001093 0.03770714
#> 4:      DS    507 1.196959 0.39776787
#> 5:      SD    607 1.164398 0.37080261
#> 6:      SO 140256 1.001093 0.03770714
#> 7:      SS    507 1.196959 0.39776787
# }
```
