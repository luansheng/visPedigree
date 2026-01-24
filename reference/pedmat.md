# Genetic Relationship Matrices and Inbreeding Coefficients

Optimized calculation of additive (A), dominance (D), epistatic (AA)
relationship matrices, their inverses, and inbreeding coefficients (f).
Uses Rcpp with Meuwissen & Luo (1992) algorithm for efficient
computation.

## Usage

``` r
pedmat(
  ped,
  method = "A",
  sparse = TRUE,
  invert_method = "auto",
  threads = 0,
  compact = FALSE
)
```

## Arguments

- ped:

  A tidied pedigree from
  [`tidyped`](https://luansheng.github.io/visPedigree/reference/tidyped.md).
  Must be a single pedigree, not a splitped object. For splitped
  results, use `pedmat(ped_split$GP1, ...)` to process individual
  groups.

- method:

  Character, one of:

  - `"A"`: Additive (numerator) relationship matrix (default)

  - `"f"`: Inbreeding coefficients (returns named vector)

  - `"Ainv"`: Inverse of A using Henderson's rules (O(n) complexity)

  - `"D"`: Dominance relationship matrix

  - `"Dinv"`: Inverse of D (requires matrix inversion)

  - `"AA"`: Additive-by-additive epistatic matrix (A \# A)

  - `"AAinv"`: Inverse of AA

- sparse:

  Logical, if `TRUE` returns sparse Matrix (recommended for large
  pedigrees). Default is `TRUE`.

- invert_method:

  Character, method for matrix inversion (Dinv/AAinv only):

  - `"auto"`: Auto-detect and use optimal method (default)

  - `"sympd"`: Force Cholesky decomposition (faster for SPD matrices)

  - `"general"`: Force general LU decomposition

- threads:

  Integer, reserved for future parallel support. Currently the C++
  implementation uses all available cores automatically.

- compact:

  Logical, if `TRUE` compacts full-sibling families by selecting one
  representative per family. This dramatically reduces matrix dimensions
  for pedigrees with large full-sib groups. See Details.

## Value

Returns a matrix or vector with S3 class `"pedmat"`.

**Object type by method:**

- `method="f"`: Named numeric vector of inbreeding coefficients

- All other methods: Sparse or dense matrix (depending on `sparse`)

**S3 Methods:**

- `print(x)`: Display matrix with metadata header

- [`summary_pedmat`](https://luansheng.github.io/visPedigree/reference/summary_pedmat.md)`(x)`:
  Detailed statistics (size, compression, mean, density)

- `dim(x)`, `length(x)`, `diag(x)`, `t(x)`: Standard operations

- `x[i, j]`: Subsetting (behaves like underlying matrix)

- `as.matrix(x)`: Convert to base matrix

**Accessing Metadata (use [`attr()`](https://rdrr.io/r/base/attr.html),
not `$`):**

- `attr(x, "ped")`: The pedigree used (or compact pedigree if
  `compact=TRUE`)

- `attr(x, "ped_compact")`: Compact pedigree (when `compact=TRUE`)

- `attr(x, "method")`: Calculation method used

- `attr(x, "call_info")`: Full calculation metadata (timing, sizes,
  etc.)

- `names(attributes(x))`: List all available attributes

**Additional attributes when `compact = TRUE`:**

- `attr(x, "compact_map")`: data.table mapping individuals to
  representatives

- `attr(x, "family_summary")`: data.table summarizing merged families

- `attr(x, "compact_stats")`: Compression statistics (ratio, n_families,
  etc.)

## Details

**API Design:**

Only a single method may be requested per call. This design prevents
accidental heavy computations. If multiple matrices are needed, call
`pedmat()` separately for each method.

**Compact Mode (`compact = TRUE`):**

Full-siblings share identical relationships with all other individuals.
Compact mode exploits this by selecting one representative per full-sib
family, dramatically reducing matrix size. For example, a pedigree with
170,000 individuals might compact to 1,800 unique relationship patterns.

Key features:

- [`query_relationship`](https://luansheng.github.io/visPedigree/reference/query_relationship.md)`(x, id1, id2)`:
  Query any individual pair, including merged siblings (automatic
  lookup)

- [`expand_pedmat`](https://luansheng.github.io/visPedigree/reference/expand_pedmat.md)`(x)`:
  Restore full matrix dimensions

- [`vismat`](https://luansheng.github.io/visPedigree/reference/vismat.md)`(x)`:
  Visualize directly (auto-detects compact)

**Performance Notes:**

- **Ainv**: O(n) complexity using Henderson's rules. Fast for any size.

- **Dinv/AAinv**: O(n³) matrix inversion. Practical limits:

  - n \< 500: ~10-20 ms

  - n = 1,000: ~40-60 ms

  - n = 2,000: ~130-150 ms

  - n \> 2,000: Consider using `compact = TRUE`

- **Memory**: Sparse matrices use ~O(nnz) memory; dense use O(n²)

## References

Meuwissen, T. H. E., & Luo, Z. (1992). Computing inbreeding coefficients
in large populations. Genetics Selection Evolution, 24(4), 305-313.

Henderson, C. R. (1976). A simple method for computing the inverse of a
numerator relationship matrix used in prediction of breeding values.
Biometrics, 32(1), 69-83.

## See also

[`tidyped`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
for preparing pedigree data,
[`query_relationship`](https://luansheng.github.io/visPedigree/reference/query_relationship.md)
for querying individual pairs,
[`expand_pedmat`](https://luansheng.github.io/visPedigree/reference/expand_pedmat.md)
for restoring full dimensions,
[`vismat`](https://luansheng.github.io/visPedigree/reference/vismat.md)
for visualization,
[`inbreed`](https://luansheng.github.io/visPedigree/reference/inbreed.md)
for simple inbreeding calculation

## Examples

``` r
# Basic usage with small pedigree
library(visPedigree)
tped <- tidyped(small_ped)

# --- Additive Relationship Matrix (default) ---
A <- pedmat(tped)
A["A", "B"]      # Relationship between A and B
#> [1] 0
diag(A)          # Diagonal = 1 + F (inbreeding)
#>        A        B        F        I       J1       J2        N        O 
#> 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 
#>        R        C        D        E        P        Q        G        H 
#> 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 
#>        K        L        M        S        T        U        V        W 
#> 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.062500 1.000000 
#>        X        Y       Z1       Z2 
#> 1.078125 1.000000 1.031250 1.031250 

# --- Inbreeding Coefficients ---
f <- pedmat(tped, method = "f")
f["Z1"]  # Inbreeding of individual Z1
#>      Z1 
#> 0.03125 

# --- Using summary_pedmat() ---
summary_pedmat(A)   # Detailed matrix statistics
#> Summary of Pedigree Matrix (A)
#> ========================================
#> Input Size:      28  individuals
#> Calculated Size: 28  individuals
#> 
#> Matrix Properties:
#> - Mean relationship:  0.352409 
#> - Density (non-zero): 54.08%
#> ========================================

# --- Accessing Metadata ---
attr(A, "ped")              # Original pedigree
#> Tidy Pedigree Object
#>        Ind   Sire    Dam   Gen    Sex IndNum SireNum DamNum Family FamilySize
#>     <char> <char> <char> <int> <char>  <int>   <int>  <int> <char>      <int>
#>  1:      A   <NA>   <NA>     1   male      1       0      0   <NA>          1
#>  2:      B   <NA>   <NA>     1 female      2       0      0   <NA>          1
#>  3:      F   <NA>   <NA>     1 female      3       0      0   <NA>          1
#>  4:      I   <NA>   <NA>     1 female      4       0      0   <NA>          1
#>  5:     J1   <NA>   <NA>     1 female      5       0      0   <NA>          1
#>  6:     J2   <NA>   <NA>     1   male      6       0      0   <NA>          1
#>  7:      N   <NA>   <NA>     1   male      7       0      0   <NA>          1
#>  8:      O   <NA>   <NA>     1 female      8       0      0   <NA>          1
#>  9:      R   <NA>   <NA>     1   male      9       0      0   <NA>          1
#> 10:      C      A      B     2 female     10       1      2    AxB          3
#> 11:      D      A      B     2   <NA>     11       1      2    AxB          3
#> 12:      E      A      B     2   male     12       1      2    AxB          3
#> 13:      P      N      O     2   <NA>     13       7      8    NxO          2
#> 14:      Q      N      O     2   male     14       7      8    NxO          2
#> 15:      G      E      F     3 female     15      12      3    ExF          2
#> 16:      H      E      F     3   male     16      12      3    ExF          2
#> 17:      K     J2      C     3 female     17       6     10   J2xC          3
#> 18:      L     J2      C     3   <NA>     18       6     10   J2xC          3
#> 19:      M     J2      C     3   male     19       6     10   J2xC          3
#> 20:      S      Q     J1     3   <NA>     20      14      5   QxJ1          2
#> 21:      T      Q     J1     3   male     21      14      5   QxJ1          2
#> 22:      U      T      K     4   male     22      21     17    TxK          1
#> 23:      V      M      G     4 female     23      19     15    MxG          1
#> 24:      W      H      I     4 female     24      16      4    HxI          1
#> 25:      X      U      V     5 female     25      22     23    UxV          1
#> 26:      Y      R      W     5   <NA>     26       9     24    RxW          1
#> 27:     Z1      N      X     6   <NA>     27       7     25    NxX          2
#> 28:     Z2      N      X     6   <NA>     28       7     25    NxX          2
#>        Ind   Sire    Dam   Gen    Sex IndNum SireNum DamNum Family FamilySize
#>     <char> <char> <char> <int> <char>  <int>   <int>  <int> <char>      <int>
attr(A, "method")           # "A"
#> [1] "A"
names(attributes(A))        # All available attributes
#>  [1] "i"         "p"         "Dim"       "Dimnames"  "x"         "uplo"     
#>  [7] "factors"   "class"     "call_info" "method"    "ped"       "pedmat_S4"

# --- Compact Mode (for large full-sib families) ---
A_compact <- pedmat(tped, method = "A", compact = TRUE)

# Query relationships (works for any individual, including merged sibs)
query_relationship(A_compact, "Z1", "Z2")  # Full-sibs Z1 and Z2
#> [1] 0.5507812

# View compression statistics
attr(A_compact, "compact_stats")
#> $n_original
#> [1] 28
#> 
#> $n_compact
#> [1] 23
#> 
#> $n_removed
#> [1] 5
#> 
#> $n_families_compacted
#> [1] 6
#> 
#> $compression_ratio
#> [1] 0.8214286
#> 
#> $memory_saved_pct
#> [1] 32.52551
#> 
#> $by_sex
#>             Sex n_original n_compact n_removed
#>          <char>      <int>     <int>     <int>
#> 1: with_parents         19        11         8
#> 2:      founder          9         9         0
#> 
#> $family_size_dist
#>    size_category n_families n_individuals_total n_individuals_removed
#>           <fctr>      <int>               <int>                 <num>
#> 1:          2-10          2                   6                     4
#> 2:             1          4                   8                     4
#> 
attr(A_compact, "family_summary")
#>    FamilyID FamilyLabel   Sire    Dam SireNum DamNum FamilySize n_compressed
#>      <char>      <char> <char> <char>   <int>  <int>      <int>        <int>
#> 1:    F0001         AxB      A      B       1      2          3            1
#> 2:    F0002         NxO      N      O       7      8          2            1
#> 3:    F0003         ExF      E      F      12      3          2            0
#> 4:    F0004        J2xC     J2      C       6     10          3            1
#> 5:    F0005        QxJ1      Q     J1      14      5          2            1
#> 6:    F0006         NxX      N      X       7     25          2            1
#>    RepInd RepIndNum   Gen
#>    <char>     <int> <int>
#> 1:      C        10     2
#> 2:      Q        12     2
#> 3:      G        13     3
#> 4:      K        15     3
#> 5:      T        17     3
#> 6:     Z1        23     6

# Expand back to full size
A_full <- expand_pedmat(A_compact)
dim(A_full)  # Original dimensions restored
#> [1] 28 28

# --- Inverse Matrices ---
Ainv <- pedmat(tped, method = "Ainv")  # Henderson's rules (fast)

# --- Dominance and Epistatic ---
D <- pedmat(tped, method = "D")
AA <- pedmat(tped, method = "AA")

# --- Visualization (requires display device) ---
if (FALSE) { # \dontrun{
vismat(A)                       # Heatmap of relationship matrix
vismat(A_compact)               # Works with compact matrices
vismat(A, grouping = "Gen")     # Group by generation
} # }
```
