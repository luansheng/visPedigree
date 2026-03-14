# Summarize Inbreeding Levels

Classifies individuals into inbreeding levels based on their inbreeding
coefficients (F) according to standard or user-defined thresholds.

## Usage

``` r
pedfclass(ped, breaks = c(0.0625, 0.125, 0.25), labels = NULL)
```

## Arguments

- ped:

  A `tidyped` object.

- breaks:

  Numeric vector of strictly increasing positive upper bounds for
  inbreeding classes. Default is `c(0.0625, 0.125, 0.25)`, corresponding
  approximately to half-sib, avuncular/grandparent, and
  full-sib/parent-offspring mating thresholds. The class `"F = 0"` is
  always kept as a fixed first level. A final open-ended class
  `"F > max(breaks)"` is always appended automatically.

- labels:

  Optional character vector of interval labels. If `NULL`, labels are
  generated automatically from `breaks`. When supplied, its length must
  equal `length(breaks)`, with each element naming the bounded interval
  `(breaks[i-1], breaks[i]]`. The open-ended tail class is always
  auto-generated and cannot be overridden.

## Value

A `data.table` with 3 columns:

- `FClass`:

  An ordered factor. By default it contains 5 levels: `"F = 0"`,
  `"0 < F <= 0.0625"`, `"0.0625 < F <= 0.125"`, `"0.125 < F <= 0.25"`,
  and `"F > 0.25"`. The number of levels equals `length(breaks) + 2`
  (the fixed zero class plus one class per bounded interval plus the
  open-ended tail).

- `Count`:

  Integer. Number of individuals in each class.

- `Percentage`:

  Numeric. Percentage of individuals in each class.

## Details

The default thresholds follow common pedigree interpretation rules:

- `F = 0.0625`: approximately the offspring of half-sib mating.

- `F = 0.125`: approximately the offspring of avuncular or
  grandparent-grandchild mating.

- `F = 0.25`: approximately the offspring of full-sib or
  parent-offspring mating.

Therefore, assigning `F = 0.25` to the class `"0.125 < F <= 0.25"` is
appropriate. If finer reporting is needed, supply custom `breaks`, for
example to separate `0.25`, `0.375`, or `0.5`.

## Examples

``` r
tp <- tidyped(simple_ped, addnum = TRUE)
pedfclass(tp)
#> Calculating inbreeding coefficients...
#> Key: <FClass>
#>                 FClass Count Percentage
#>                  <ord> <int>      <num>
#> 1:               F = 0    58  98.305085
#> 2:     0 < F <= 0.0625     1   1.694915
#> 3: 0.0625 < F <= 0.125     0   0.000000
#> 4:   0.125 < F <= 0.25     0   0.000000
#> 5:            F > 0.25     0   0.000000

# Finer custom classes (4 breaks, labels auto-generated)
pedfclass(tp, breaks = c(0.03125, 0.0625, 0.125, 0.25))
#> Calculating inbreeding coefficients...
#> Key: <FClass>
#>                   FClass Count Percentage
#>                    <ord> <int>      <num>
#> 1:                 F = 0    58  98.305085
#> 2:      0 < F <= 0.03125     1   1.694915
#> 3: 0.03125 < F <= 0.0625     0   0.000000
#> 4:   0.0625 < F <= 0.125     0   0.000000
#> 5:     0.125 < F <= 0.25     0   0.000000
#> 6:              F > 0.25     0   0.000000

# Custom labels aligned to breaks (3 labels for 3 breaks; tail is auto)
pedfclass(tp, labels = c("Low", "Moderate", "High"))
#> Calculating inbreeding coefficients...
#> Key: <FClass>
#>      FClass Count Percentage
#>       <ord> <int>      <num>
#> 1:    F = 0    58  98.305085
#> 2:      Low     1   1.694915
#> 3: Moderate     0   0.000000
#> 4:     High     0   0.000000
#> 5: F > 0.25     0   0.000000

# \donttest{
tp_inbred <- tidyped(inbred_ped, addnum = TRUE)
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
# }
```
