# Calculate Equi-Generate Coefficient

Estimates the number of distinct ancestral generations using the
Equi-Generate Coefficient (ECG). The ECG is calculated as 1/2 of the sum
of the parents' ECG values plus 1.

## Usage

``` r
pedecg(ped)
```

## Arguments

- ped:

  A `tidyped` object.

## Value

A `data.table` with columns:

- `Ind`: Individual ID.

- `ECG`: Equi-Generate Coefficient.

- `FullGen`: Fully traced generations (depth of complete ancestry).

- `MaxGen`: Maximum ancestral generations (depth of deepest ancestor).

## References

Boichard, D., Maignel, L., & Verrier, E. (1997). The value of using
probabilities of gene origin to measure genetic variability in a
population. Genetics Selection Evolution, 29(1), 5.

## Examples

``` r
tp <- tidyped(simple_ped)
ecg <- pedecg(tp)

# ECG combines pedigree depth and completeness
head(ecg)
#>       Ind   ECG FullGen MaxGen
#>    <char> <num>   <num>  <num>
#> 1: J0C032     0       0      0
#> 2: J0C185     0       0      0
#> 3: J0C231     0       0      0
#> 4: J0C317     0       0      0
#> 5: J0C355     0       0      0
#> 6: J0C450     0       0      0

# Individuals with deeper and more complete ancestry have larger ECG values
ecg[order(-ECG)][1:5]
#>       Ind     ECG FullGen MaxGen
#>    <char>   <num>   <num>  <num>
#> 1: J5X804 4.46875       2      5
#> 2: J4Y326 3.93750       3      4
#> 3: J3Y620 3.00000       3      3
#> 4: J4E185 3.00000       1      4
#> 5: J3Y771 2.87500       2      3
```
