# Calculate Generation Intervals

Computes the generation intervals for the four gametic pathways: Sire to
Son (SS), Sire to Daughter (SD), Dam to Son (DS), and Dam to Daughter
(DD). The generation interval is defined as the age of the parents at
the birth of their offspring.

## Usage

``` r
pedgenint(
  ped,
  timevar = NULL,
  unit = c("year", "month", "day", "hour"),
  format = NULL,
  cycle = NULL,
  by = NULL
)
```

## Arguments

- ped:

  A `tidyped` object.

- timevar:

  Character. The name of the column containing the **birth date** (or
  hatch date) of each individual. The column must be one of:

  - `Date` or `POSIXct` (recommended).

  - A date string parseable by `as.POSIXct` (e.g., `"2020-03-15"`). Use
    `format` for non-ISO strings.

  - A numeric year (e.g., `2020`). Automatically converted to `Date`
    (`"YYYY-07-01"`) with a message. For finer precision, convert to
    `Date` beforehand.

  If `NULL`, auto-detects columns named `"BirthYear"`, `"Year"`,
  `"BirthDate"`, or `"Date"`.

- unit:

  Character. Output time unit for the interval: `"year"` (default),
  `"month"`, `"day"`, or `"hour"`.

- format:

  Character. Optional format string for parsing `timevar` when it
  contains non-standard date strings (e.g., `"%d/%m/%Y"` for
  `"15/03/2020"`).

- cycle:

  Numeric. Optional target (designed) length of one generation cycle
  expressed in `unit`s. When provided, an additional column `GenEquiv`
  is appended to the result, defined as: \$\$GenEquiv_i =
  \frac{\bar{L}\_i}{L\_{cycle}}\$\$ where \\\bar{L}\_i\\ is the observed
  mean interval for pathway \\i\\ and \\L\_{cycle}\\ is `cycle`. A value
  \> 1 means the observed interval exceeds the target cycle (lower
  breeding efficiency). Example: for Pacific white shrimp with a 180-day
  target cycle, set `unit = "day", cycle = 180`.

- by:

  Character. Optional grouping column (e.g., `"Breed"`, `"Farm"`). If
  provided, intervals are calculated within each group.

## Value

A `data.table` with columns:

- `Group`: Grouping level (if `by` is used).

- `Pathway`: One of `"SS"`, `"SD"`, `"DS"`, `"DD"`, `"SO"`, `"DO"`, or
  `"Average"`. SS/SD/DS/DD require offspring sex; SO (Sire-Offspring)
  and DO (Dam-Offspring) are computed from all parent-offspring pairs
  regardless of offspring sex.

- `N`: Number of parent-offspring pairs used.

- `Mean`: Average generation interval in `unit`.

- `SD`: Standard deviation of the interval.

- `GenEquiv`: (Optional) Generation equivalents based on `cycle`.

## Details

Parent-offspring pairs with zero or negative intervals are excluded from
the calculation because they typically indicate data entry errors or
insufficient time resolution. If many zero intervals are expected (e.g.,
when using `unit = "year"` with annual spawners), consider using a finer
time unit such as `"month"` or `"day"`.

Numeric year columns (e.g., `2020`) are automatically converted to
`Date` by appending `"-07-01"` (mid-year) as a reasonable default. For
more precise results, convert to `Date` before calling this function.

## Examples

``` r
# \donttest{
# ---- Basic usage with package dataset (numeric Year auto-converted) ----
tped <- tidyped(big_family_size_ped)
gi <- pedgenint(tped, timevar = "Year")
#> Numeric time column detected. Converting to Date (YYYY-07-01). For finer precision, convert to Date beforehand.
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
gi
#>    Pathway      N     Mean         SD
#>     <char>  <int>    <num>      <num>
#> 1: Average 280512 1.001093 0.03770707
#> 2:      DD    607 1.164398 0.37080261
#> 3:      DO 140256 1.001093 0.03770714
#> 4:      DS    507 1.196959 0.39776787
#> 5:      SD    607 1.164398 0.37080261
#> 6:      SO 140256 1.001093 0.03770714
#> 7:      SS    507 1.196959 0.39776787

# ---- Generation equivalents with cycle ----
gi2 <- pedgenint(tped, timevar = "Year", cycle = 2)
#> Numeric time column detected. Converting to Date (YYYY-07-01). For finer precision, convert to Date beforehand.
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
#> Warning: Subsetting removed parent records. Result is a plain data.table, not a tidyped.
#> Use tidyped(tp, cand = ids, trace = "up") to extract a valid sub-pedigree.
gi2
#>    Pathway      N     Mean         SD  GenEquiv
#>     <char>  <int>    <num>      <num>     <num>
#> 1: Average 280512 1.001093 0.03770707 0.5005467
#> 2:      DD    607 1.164398 0.37080261 0.5821992
#> 3:      DO 140256 1.001093 0.03770714 0.5005467
#> 4:      DS    507 1.196959 0.39776787 0.5984796
#> 5:      SD    607 1.164398 0.37080261 0.5821992
#> 6:      SO 140256 1.001093 0.03770714 0.5005467
#> 7:      SS    507 1.196959 0.39776787 0.5984796
# }
```
