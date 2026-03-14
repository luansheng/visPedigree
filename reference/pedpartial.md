# Calculate Partial Inbreeding

Decomposes individuals' inbreeding coefficients into marginal
contributions from specific ancestors. This allows identifying which
ancestors or lineages are responsible for the observed inbreeding.

## Usage

``` r
pedpartial(ped, ancestors = NULL, top = 20)
```

## Arguments

- ped:

  A `tidyped` object.

- ancestors:

  Character vector. IDs of ancestors to calculate partial inbreeding
  for. If NULL, the top ancestors by marginal contribution are used.

- top:

  Integer. Number of top ancestors to include if `ancestors` is NULL.

## Value

A `data.table` with the first column as `Ind` and subsequent columns
representing the partial inbreeding (\$pF\$) from each ancestor.

## Details

The sum of all partial inbreeding coefficients for an individual
(including contributions from founders) equals \$1 + f_i\$, where
\$f_i\$ is the total inbreeding coefficient. This function specifically
isolates the terms in the Meuwissen & Luo (1992) decomposition that
correspond to the selected ancestors.

## References

Lacey, R. C. (1996). A formula for determining the partial inbreeding
coefficient, \$F_ij\$. Journal of Heredity, 87(4), 337-339.

Meuwissen, T. H., & Luo, Z. (1992). Computing inbreeding coefficients in
large populations. Genetics Selection Evolution, 24(4), 305-313.

## Examples

``` r
# \donttest{
library(data.table)
tp <- tidyped(inbred_ped)
# Calculate partial inbreeding originating from specific ancestors
target_ancestors <- inbred_ped[is.na(Sire) & is.na(Dam), Ind]
pF <- pedpartial(tp, ancestors = target_ancestors)
#> Calculating partial inbreeding for 2 ancestors...
print(tail(pF))
#>       Ind       A       B
#>    <char>   <num>   <num>
#> 1:      B 0.00000 0.00000
#> 2:      C 0.00000 0.00000
#> 3:      D 0.00000 0.00000
#> 4:      E 0.12500 0.12500
#> 5:      F 0.00000 0.25000
#> 6:      G 0.15625 0.28125
# }
```
