# Split Pedigree into Disconnected Groups

Detects and splits a tidyped object into disconnected groups (connected
components). Uses igraph to efficiently find groups of individuals that
have no genetic relationships with each other. Isolated individuals (Gen
= 0, those with no parents and no offspring) are excluded from group
splitting and stored separately.

## Usage

``` r
splitped(ped)
```

## Arguments

- ped:

  A tidyped object created by
  [`tidyped`](https://luansheng.github.io/visPedigree/reference/tidyped.md).

## Value

A list of class "splitped" containing:

- GP1, GP2, ...:

  tidyped objects for each disconnected group (with at least 2
  individuals), with renumbered IndNum, SireNum, DamNum

The returned object has the following attributes:

- n_groups:

  Number of disconnected groups found (excluding isolated individuals)

- sizes:

  Named vector of group sizes

- total:

  Total number of individuals in groups (excluding isolated)

- isolated:

  Character vector of isolated individual IDs (Gen = 0)

- n_isolated:

  Number of isolated individuals

## Details

This function identifies connected components in the pedigree graph
where edges represent parent-offspring relationships. Two individuals
are in the same group if they share any ancestry (direct or indirect).

Isolated individuals (Gen = 0 in tidyped output) are those who:

- Have no known parents (Sire and Dam are both NA)

- Are not parents of any other individual in the pedigree

These isolated individuals are excluded from splitting and stored in the
`isolated` attribute. Each resulting group contains at least 2
individuals (at least one parent-offspring relationship).

The function always returns a list, even if there is only one group
(i.e., the pedigree is fully connected). Groups are sorted by size in
descending order.

Each group in the result is a valid tidyped object with:

- Renumbered IndNum (1 to n for each group)

- Updated SireNum and DamNum referencing the new IndNum

- Recalculated Gen (generation) based on the group's structure

## See also

[`tidyped`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
for pedigree tidying

## Examples

``` r
# Load example data
library(visPedigree)
data(small_ped)

# First tidy the pedigree
tped <- tidyped(small_ped)

# Split into groups
result <- splitped(tped)
print(result)
#> Pedigree Split Result
#> ======================
#> Total individuals in groups: 28 
#> Isolated individuals (Gen=0): 0 
#> Number of groups:  1 
#> 
#> Group sizes:
#>   GP1: 28 individuals (100.0%)

# Access individual groups (each is a tidyped object)
result$GP1
#> Tidy Pedigree Object
#> Indices: <Sire>, <Dam>, <Sex>
#>        Ind   Sire    Dam    Sex   Gen IndNum SireNum DamNum Family FamilySize
#>     <char> <char> <char> <char> <int>  <int>   <int>  <int> <char>      <int>
#>  1:      A   <NA>   <NA>   male     1      1       0      0   <NA>          1
#>  2:      B   <NA>   <NA> female     1      2       0      0   <NA>          1
#>  3:      F   <NA>   <NA> female     1      3       0      0   <NA>          1
#>  4:      I   <NA>   <NA> female     1      4       0      0   <NA>          1
#>  5:     J1   <NA>   <NA> female     1      5       0      0   <NA>          1
#>  6:     J2   <NA>   <NA>   male     1      6       0      0   <NA>          1
#>  7:      N   <NA>   <NA>   male     1      7       0      0   <NA>          1
#>  8:      O   <NA>   <NA> female     1      8       0      0   <NA>          1
#>  9:      R   <NA>   <NA>   male     1      9       0      0   <NA>          1
#> 10:      C      A      B female     2     10       1      2    AxB          3
#> 11:      D      A      B   <NA>     2     11       1      2    AxB          3
#> 12:      E      A      B   male     2     12       1      2    AxB          3
#> 13:      P      N      O   <NA>     2     13       7      8    NxO          2
#> 14:      Q      N      O   male     2     14       7      8    NxO          2
#> 15:      G      E      F female     3     15      12      3    ExF          2
#> 16:      H      E      F   male     3     16      12      3    ExF          2
#> 17:      K     J2      C female     3     17       6     10   J2xC          3
#> 18:      L     J2      C   <NA>     3     18       6     10   J2xC          3
#> 19:      M     J2      C   male     3     19       6     10   J2xC          3
#> 20:      S      Q     J1   <NA>     3     20      14      5   QxJ1          2
#> 21:      T      Q     J1   male     3     21      14      5   QxJ1          2
#> 22:      U      T      K   male     4     22      21     17    TxK          1
#> 23:      V      M      G female     4     23      19     15    MxG          1
#> 24:      W      H      I female     4     24      16      4    HxI          1
#> 25:      X      U      V female     5     25      22     23    UxV          1
#> 26:      Y      R      W   <NA>     5     26       9     24    RxW          1
#> 27:     Z1      N      X   <NA>     6     27       7     25    NxX          2
#> 28:     Z2      N      X   <NA>     6     28       7     25    NxX          2
#>        Ind   Sire    Dam    Sex   Gen IndNum SireNum DamNum Family FamilySize
#>     <char> <char> <char> <char> <int>  <int>   <int>  <int> <char>      <int>

# Check isolated individuals
attr(result, "isolated")
#> character(0)
```
