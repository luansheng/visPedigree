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

  Integer. Number of OpenMP threads to use. Use 0 to keep the
  system/default setting. Currently, multi-threading is explicitly
  implemented for:

  - `"D"`: Dominance relationship matrix (significant speedup).

  - `"Ainv"`: Inverse of A (only for large pedigrees, n \>= 5000).

  For `"Dinv"`, `"AA"`, and `"AAinv"`, parallelism depends on the linked
  BLAS/LAPACK library (e.g., OpenBLAS, MKL, Accelerate) and is not
  controlled by this parameter. Methods `"A"` and `"f"` are
  single-threaded.

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

- `dim(x)`, `length(x)`, `Matrix::diag(x)`, `t(x)`: Standard operations

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
  Visualize directly (auto-expands compact)

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
[`pedprod`](https://luansheng.github.io/visPedigree/reference/pedprod.md)
for matrix-free products with A or Ainv,
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
Matrix::diag(A)  # Diagonal = 1 + F (inbreeding)
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
#> - Mean off-diagonal relationship:  0.136088 
#> - Density (non-zero): 54.08%
#> ========================================

# --- Accessing Metadata ---
attr(A, "ped")              # Original pedigree
#> Tidy Pedigree Object
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
attr(A, "method")           # "A"
#> [1] "A"
names(attributes(A))        # All available attributes
#> [1] "Dim"       "Dimnames"  "x"         "factors"   "class"     "call_info"
#> [7] "method"    "ped"       "pedmat_S4"

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
#> [1] 27
#> 
#> $n_removed
#> [1] 1
#> 
#> $n_families_compacted
#> [1] 5
#> 
#> $compression_ratio
#> [1] 0.9642857
#> 
#> $memory_saved_pct
#> [1] 7.015306
#> 
#> $by_sex
#>             Sex n_original n_compact n_removed
#>          <char>      <int>     <int>     <int>
#> 1: with_parents         19        18         1
#> 2:      founder          9         9         0
#> 
#> $family_size_dist
#>    SizeCategory n_families n_individuals_total n_individuals_removed
#>          <fctr>      <int>               <int>                 <num>
#> 1:            1          5                   6                     1
#> 
attr(A_compact, "family_summary")
#>    FamilyID FamilyLabel   Sire    Dam SireNum DamNum FamilySize NCompressed
#>      <char>      <char> <char> <char>   <int>  <int>      <int>       <int>
#> 1:    F0001         AxB      A      B       1      2          1           0
#> 2:    F0002         NxO      N      O       7      8          1           0
#> 3:    F0003        J2xC     J2      C       6     10          1           0
#> 4:    F0004        QxJ1      Q     J1      14      5          1           0
#> 5:    F0005         NxX      N      X       7     25          2           1
#>    RepInd RepIndNum   Gen
#>    <char>     <int> <int>
#> 1:      D        11     2
#> 2:      P        13     2
#> 3:      L        18     3
#> 4:      S        20     3
#> 5:     Z1        27     6

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
vismat(A, by = "Gen")     # Group by generation
} # }
```
