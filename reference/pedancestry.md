# Calculate Ancestry Proportions

Estimates the proportion of genes for each individual that originates
from specific founder groups (e.g., breeds, source populations).

## Usage

``` r
pedancestry(ped, foundervar, target_labels = NULL)
```

## Arguments

- ped:

  A `tidyped` object.

- foundervar:

  Character. The name of the column containing founder-group labels
  (e.g., "Breed", "Origin").

- target_labels:

  Character vector. Specific founder-group labels to track. If NULL, all
  unique labels in `foundervar` among founders are used.

## Value

A `data.table` with columns:

- `Ind`: Individual ID.

- One column per tracked label (named after each unique value in
  `foundervar` among founders, or as specified by `target_labels`). Each
  value gives the proportion of genes (0–1) originating from that
  founder group. Row sums across all label columns equal 1.

## Examples

``` r
# \donttest{
library(data.table)
#> 
#> Attaching package: ‘data.table’
#> The following object is masked from ‘package:base’:
#> 
#>     %notin%
# Create dummy labels for founders
tp <- tidyped(small_ped)
tp_dated <- copy(tp)
founders <- tp_dated[is.na(Sire) & is.na(Dam), Ind]
# Assign 'LineA' and 'LineB'
tp_dated[Ind %in% founders[1:(length(founders)/2)], Origin := "LineA"]
#> Tidy Pedigree Object
#> Index: <Ind>
#>        Ind   Sire    Dam Family FamilySize   Gen    Sex IndNum SireNum DamNum
#>     <char> <char> <char> <char>      <int> <int> <char>  <int>   <int>  <int>
#>  1:      A   <NA>   <NA>   <NA>          1     1   male      1       0      0
#>  2:      B   <NA>   <NA>   <NA>          1     1 female      2       0      0
#>  3:      F   <NA>   <NA>   <NA>          1     1 female      3       0      0
#>  4:      I   <NA>   <NA>   <NA>          1     1 female      4       0      0
#>  5:     J1   <NA>   <NA>   <NA>          1     1 female      5       0      0
#>  6:     J2   <NA>   <NA>   <NA>          1     1   male      6       0      0
#>  7:      N   <NA>   <NA>   <NA>          1     1   male      7       0      0
#>  8:      O   <NA>   <NA>   <NA>          1     1 female      8       0      0
#>  9:      R   <NA>   <NA>   <NA>          1     1   male      9       0      0
#> 10:      C      A      B    AxB          3     2 female     10       1      2
#> 11:      D      A      B    AxB          3     2   <NA>     11       1      2
#> 12:      E      A      B    AxB          3     2   male     12       1      2
#> 13:      P      N      O    NxO          2     2   <NA>     13       7      8
#> 14:      Q      N      O    NxO          2     2   male     14       7      8
#> 15:      G      E      F    ExF          2     3 female     15      12      3
#> 16:      H      E      F    ExF          2     3   male     16      12      3
#> 17:      K     J2      C   J2xC          3     3 female     17       6     10
#> 18:      L     J2      C   J2xC          3     3   <NA>     18       6     10
#> 19:      M     J2      C   J2xC          3     3   male     19       6     10
#> 20:      S      Q     J1   QxJ1          2     3   <NA>     20      14      5
#> 21:      T      Q     J1   QxJ1          2     3   male     21      14      5
#> 22:      U      T      K    TxK          1     4   male     22      21     17
#> 23:      V      M      G    MxG          1     4 female     23      19     15
#> 24:      W      H      I    HxI          1     4 female     24      16      4
#> 25:      X      U      V    UxV          1     5 female     25      22     23
#> 26:      Y      R      W    RxW          1     5   <NA>     26       9     24
#> 27:     Z1      N      X    NxX          2     6   <NA>     27       7     25
#> 28:     Z2      N      X    NxX          2     6   <NA>     28       7     25
#>        Ind   Sire    Dam Family FamilySize   Gen    Sex IndNum SireNum DamNum
#>     <char> <char> <char> <char>      <int> <int> <char>  <int>   <int>  <int>
#>     Origin
#>     <char>
#>  1:  LineA
#>  2:  LineA
#>  3:  LineA
#>  4:  LineA
#>  5:   <NA>
#>  6:   <NA>
#>  7:   <NA>
#>  8:   <NA>
#>  9:   <NA>
#> 10:   <NA>
#> 11:   <NA>
#> 12:   <NA>
#> 13:   <NA>
#> 14:   <NA>
#> 15:   <NA>
#> 16:   <NA>
#> 17:   <NA>
#> 18:   <NA>
#> 19:   <NA>
#> 20:   <NA>
#> 21:   <NA>
#> 22:   <NA>
#> 23:   <NA>
#> 24:   <NA>
#> 25:   <NA>
#> 26:   <NA>
#> 27:   <NA>
#> 28:   <NA>
#>     Origin
#>     <char>
tp_dated[is.na(Origin), Origin := "LineB"]
#> Tidy Pedigree Object
#> Index: <Ind>
#>        Ind   Sire    Dam Family FamilySize   Gen    Sex IndNum SireNum DamNum
#>     <char> <char> <char> <char>      <int> <int> <char>  <int>   <int>  <int>
#>  1:      A   <NA>   <NA>   <NA>          1     1   male      1       0      0
#>  2:      B   <NA>   <NA>   <NA>          1     1 female      2       0      0
#>  3:      F   <NA>   <NA>   <NA>          1     1 female      3       0      0
#>  4:      I   <NA>   <NA>   <NA>          1     1 female      4       0      0
#>  5:     J1   <NA>   <NA>   <NA>          1     1 female      5       0      0
#>  6:     J2   <NA>   <NA>   <NA>          1     1   male      6       0      0
#>  7:      N   <NA>   <NA>   <NA>          1     1   male      7       0      0
#>  8:      O   <NA>   <NA>   <NA>          1     1 female      8       0      0
#>  9:      R   <NA>   <NA>   <NA>          1     1   male      9       0      0
#> 10:      C      A      B    AxB          3     2 female     10       1      2
#> 11:      D      A      B    AxB          3     2   <NA>     11       1      2
#> 12:      E      A      B    AxB          3     2   male     12       1      2
#> 13:      P      N      O    NxO          2     2   <NA>     13       7      8
#> 14:      Q      N      O    NxO          2     2   male     14       7      8
#> 15:      G      E      F    ExF          2     3 female     15      12      3
#> 16:      H      E      F    ExF          2     3   male     16      12      3
#> 17:      K     J2      C   J2xC          3     3 female     17       6     10
#> 18:      L     J2      C   J2xC          3     3   <NA>     18       6     10
#> 19:      M     J2      C   J2xC          3     3   male     19       6     10
#> 20:      S      Q     J1   QxJ1          2     3   <NA>     20      14      5
#> 21:      T      Q     J1   QxJ1          2     3   male     21      14      5
#> 22:      U      T      K    TxK          1     4   male     22      21     17
#> 23:      V      M      G    MxG          1     4 female     23      19     15
#> 24:      W      H      I    HxI          1     4 female     24      16      4
#> 25:      X      U      V    UxV          1     5 female     25      22     23
#> 26:      Y      R      W    RxW          1     5   <NA>     26       9     24
#> 27:     Z1      N      X    NxX          2     6   <NA>     27       7     25
#> 28:     Z2      N      X    NxX          2     6   <NA>     28       7     25
#>        Ind   Sire    Dam Family FamilySize   Gen    Sex IndNum SireNum DamNum
#>     <char> <char> <char> <char>      <int> <int> <char>  <int>   <int>  <int>
#>     Origin
#>     <char>
#>  1:  LineA
#>  2:  LineA
#>  3:  LineA
#>  4:  LineA
#>  5:  LineB
#>  6:  LineB
#>  7:  LineB
#>  8:  LineB
#>  9:  LineB
#> 10:  LineB
#> 11:  LineB
#> 12:  LineB
#> 13:  LineB
#> 14:  LineB
#> 15:  LineB
#> 16:  LineB
#> 17:  LineB
#> 18:  LineB
#> 19:  LineB
#> 20:  LineB
#> 21:  LineB
#> 22:  LineB
#> 23:  LineB
#> 24:  LineB
#> 25:  LineB
#> 26:  LineB
#> 27:  LineB
#> 28:  LineB
#>     Origin
#>     <char>

# Calculate ancestry proportions for all individuals
anc <- pedancestry(tp_dated, foundervar = "Origin")
print(tail(anc))
#>       Ind LineA LineB
#>    <char> <num> <num>
#> 1:      V  0.75  0.25
#> 2:      W  1.00  0.00
#> 3:      X  0.50  0.50
#> 4:      Y  0.50  0.50
#> 5:     Z1  0.25  0.75
#> 6:     Z2  0.25  0.75
# }
```
