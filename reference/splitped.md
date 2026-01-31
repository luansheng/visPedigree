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
#>        Ind   Sire    Dam    Sex Family FamilySize   Gen IndNum SireNum DamNum
#>     <char> <char> <char> <char> <char>      <int> <int>  <int>   <int>  <int>
#>  1:      A   <NA>   <NA>   male   <NA>          1     1      1       0      0
#>  2:      B   <NA>   <NA> female   <NA>          1     1      2       0      0
#>  3:      F   <NA>   <NA> female   <NA>          1     1      3       0      0
#>  4:      I   <NA>   <NA> female   <NA>          1     1      4       0      0
#>  5:     J1   <NA>   <NA> female   <NA>          1     1      5       0      0
#>  6:     J2   <NA>   <NA>   male   <NA>          1     1      6       0      0
#>  7:      N   <NA>   <NA>   male   <NA>          1     1      7       0      0
#>  8:      O   <NA>   <NA> female   <NA>          1     1      8       0      0
#>  9:      R   <NA>   <NA>   male   <NA>          1     1      9       0      0
#> 10:      C      A      B female    AxB          3     2     10       1      2
#> 11:      D      A      B   <NA>    AxB          3     2     11       1      2
#> 12:      E      A      B   male    AxB          3     2     12       1      2
#> 13:      P      N      O   <NA>    NxO          2     2     13       7      8
#> 14:      Q      N      O   male    NxO          2     2     14       7      8
#> 15:      G      E      F female    ExF          2     3     15      12      3
#> 16:      H      E      F   male    ExF          2     3     16      12      3
#> 17:      K     J2      C female   J2xC          3     3     17       6     10
#> 18:      L     J2      C   <NA>   J2xC          3     3     18       6     10
#> 19:      M     J2      C   male   J2xC          3     3     19       6     10
#> 20:      S      Q     J1   <NA>   QxJ1          2     3     20      14      5
#> 21:      T      Q     J1   male   QxJ1          2     3     21      14      5
#> 22:      U      T      K   male    TxK          1     4     22      21     17
#> 23:      V      M      G female    MxG          1     4     23      19     15
#> 24:      W      H      I female    HxI          1     4     24      16      4
#> 25:      X      U      V female    UxV          1     5     25      22     23
#> 26:      Y      R      W   <NA>    RxW          1     5     26       9     24
#> 27:     Z1      N      X   <NA>    NxX          2     6     27       7     25
#> 28:     Z2      N      X   <NA>    NxX          2     6     28       7     25
#>        Ind   Sire    Dam    Sex Family FamilySize   Gen IndNum SireNum DamNum
#>     <char> <char> <char> <char> <char>      <int> <int>  <int>   <int>  <int>

# Check isolated individuals
attr(result, "isolated")
#> character(0)
```
