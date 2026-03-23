# 5. Calculation and visualization of relationship matrix

1.  [Calculating Relationship Matrices with pedmat()](#id_1)  
    1.1 [Supported Methods](#id_1-1)  
    1.2 [Basic Usage](#id_1-2)  
    1.3 [Sparse Matrix Representation](#id_1-3)  
2.  [Inspecting the Matrix](#id_2)  
    2.1 [Summary Statistics](#id_2-1)  
    2.2 [Querying Specific Relationships](#id_2-2)  
3.  [Compact Mode for Large Pedigrees](#id_3)  
    3.1 [Using compact = TRUE](#id_3-1)  
    3.2 [Expanding and Querying Compacted Matrices](#id_3-2)  
    3.3 [When to Use Compact Mode](#id_3-3)  
4.  [Visualizing Relationship Matrices with vismat()](#id_4)  
    4.1 [Relationship Heatmaps](#id_4-1)  
    4.2 [Inbreeding and Kinship Histograms](#id_4-2)  
5.  [Performance Considerations](#id_5)

Relationship matrices are fundamental tools in quantitative genetics and
animal breeding. They quantify the genetic similarity between
individuals due to shared ancestry, which is essential for estimating
breeding values (BLUP) and managing genetic diversity. The `visPedigree`
package provides efficient tools for calculating various relationship
matrices and visualizing them through heatmaps and histograms.

## 1. Calculating Relationship Matrices with `pedmat()`

The
[`pedmat()`](https://luansheng.github.io/visPedigree/reference/pedmat.md)
function is the primary tool for calculating relationship matrices. It
supports both additive and dominance relationship matrices, as well as
their inverses.

### 1.1 Supported Methods

The `method` parameter in
[`pedmat()`](https://luansheng.github.io/visPedigree/reference/pedmat.md)
determines the type of matrix to calculate:

- **“A”**: Additive relationship matrix (Numerator Relationship Matrix).
- **“Ainv”**: Inverse of the additive relationship matrix.
- **“D”**: Dominance relationship matrix.
- **“Dinv”**: Inverse of the dominance relationship matrix.
- **“AA”**: Additive-by-additive (epistatic) relationship matrix.
- **“AAinv”**: Inverse of the epistatic relationship matrix.
- **“f”**: Inbreeding coefficients vector (uses the same optimized
  engine as `tidyped(..., inbreed = TRUE)`).

### 1.2 Basic Usage

Most calculations require a pedigree tidied by
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md).

``` r
# Load example pedigree and tidy it
data(small_ped)
tped <- tidyped(small_ped)

# Calculate Additive Relationship Matrix (A)
mat_A <- pedmat(tped, method = "A")

# Calculate Dominance Relationship Matrix (D)
mat_D <- pedmat(tped, method = "D")

# Calculate inbreeding coefficients (f)
vec_f <- pedmat(tped, method = "f")
```

### 1.3 Sparse Matrix Representation

By default,
[`pedmat()`](https://luansheng.github.io/visPedigree/reference/pedmat.md)
returns a sparse matrix (class `dsCMatrix` from the `Matrix` package)
for relationship matrices. This is highly memory-efficient for large
pedigrees where many individuals are unrelated.

``` r
class(mat_A)
#> [1] "dsCMatrix"
#> attr(,"package")
#> [1] "Matrix"
```

## 2. Inspecting the Matrix

### 2.1 Summary Statistics

Use the [`summary()`](https://rdrr.io/r/base/summary.html) method to get
an overview of the calculated matrix, including size, density, and
average relationship.

``` r
tail(summary(mat_A),10)
#> 28 x 28 sparse Matrix of class "dsCMatrix", with 226 entries
#>      i  j         x
#> 217 19 28 0.2031250
#> 218 20 28 0.1875000
#> 219 21 28 0.2500000
#> 220 22 28 0.3515625
#> 221 23 28 0.3046875
#> 222 24 28 0.0468750
#> 223 25 28 0.5703125
#> 224 26 28 0.0234375
#> 225 27 28 0.5507812
#> 226 28 28 1.0312500
```

### 2.2 Querying Specific Relationships

Instead of manually indexing the matrix, you can use
[`query_relationship()`](https://luansheng.github.io/visPedigree/reference/query_relationship.md)
to retrieve coefficients by individual IDs.

``` r
# Query relationship between Z1 and Z2
query_relationship(mat_A, "Z1", "Z2")
#> [1] 0.5507812

# Query multiple pairs
query_relationship(mat_A, c("Z1", "A"), c("Z2", "B"))
#> 2 x 2 sparse Matrix of class "dgCMatrix"
#>           Z2       B
#> Z1 0.5507812 0.09375
#> A  0.0937500 .
```

## 3. Compact Mode for Large Pedigrees

For large pedigrees with many full-sibling families (common in aquatic
breeding populations),
[`pedmat()`](https://luansheng.github.io/visPedigree/reference/pedmat.md)
can merge full siblings into representative nodes to save memory and
time.

### 3.1 Using `compact = TRUE`

When `compact = TRUE`, the matrix is calculated for unique
representative individuals from each full-sib family.

``` r
# Calculate compacted A matrix
mat_compact <- pedmat(tped, method = "A", compact = TRUE)

# The result is a 'pedmat' object containing the compacted matrix
print(mat_compact[11:20,11:20])
#> 10 x 10 sparse Matrix of class "dsCMatrix"
#>   [[ suppressing 10 column names 'D', 'E', 'P' ... ]]
#>                                                        
#> D 1.00 0.50 .    .   0.250 0.250 0.250 0.250 0.250 .   
#> E 0.50 1.00 .    .   0.500 0.500 0.250 0.250 0.250 .   
#> P .    .    1.00 0.5 .     .     .     .     .     0.25
#> Q .    .    0.50 1.0 .     .     .     .     .     0.50
#> G 0.25 0.50 .    .   1.000 0.500 0.125 0.125 0.125 .   
#> H 0.25 0.50 .    .   0.500 1.000 0.125 0.125 0.125 .   
#> K 0.25 0.25 .    .   0.125 0.125 1.000 0.500 0.500 .   
#> L 0.25 0.25 .    .   0.125 0.125 0.500 1.000 0.500 .   
#> M 0.25 0.25 .    .   0.125 0.125 0.500 0.500 1.000 .   
#> S .    .    0.25 0.5 .     .     .     .     .     1.00
```

### 3.2 Expanding and Querying Compacted Matrices

If you need the full matrix after a compact calculation, use
[`expand_pedmat()`](https://luansheng.github.io/visPedigree/reference/expand_pedmat.md).
For retrieving specific values,
[`query_relationship()`](https://luansheng.github.io/visPedigree/reference/query_relationship.md)
handles both standard and compact objects transparently.

``` r
# Expand to full 28x28 matrix
mat_full <- expand_pedmat(mat_compact)
dim(mat_full)
#> [1] 28 28

# Query still works the same way
query_relationship(mat_compact, "Z1", "Z2")
#> [1] 0.5507812
```

### 3.3 When to Use Compact Mode

Compact mode is highly recommended for:

- **Large Pedigrees**: More than 5,000 individuals with substantial
  full-sibling groups.
- **High-fecundity species**: Such as aquatic animals or plants, where
  families often have hundreds or thousands of offspring.
- **Memory-limited environments**: When the full matrix exceeds
  available RAM.

| Pedigree Size | Full-Sib Proportion | Recommended Mode   |
|:--------------|:--------------------|:-------------------|
| \< 1,000      | Any                 | Standard           |
| \> 5,000      | \< 20%              | Standard / Compact |
| \> 5,000      | \> 20%              | **Compact**        |

## 4. Visualizing Relationship Matrices with `vismat()`

Visualization helps in understanding population structure, detecting
family clusters, and checking the distribution of genetic relationships.

### 4.1 Relationship Heatmaps

The `"heatmap"` type (default) uses a Nature Genetics style color
palette (White–Orange–Red) to display relationships. Rows and columns
are reordered by hierarchical clustering (Ward.D2) by default, bringing
closely related individuals into contiguous blocks — full-sibs cluster
tightly because they share nearly identical relationship profiles with
the rest of the population.

``` r
# Heatmap of the A matrix (with default clustering reorder)
vismat(mat_A, labelcex = 0.5)
```

![](relationship-matrix_files/figure-html/heatmap-1.png)

#### Compact Matrix — Direct Visualization

A compact `pedmat` object can be passed directly to
[`vismat()`](https://luansheng.github.io/visPedigree/reference/vismat.md).
It is automatically expanded to full dimensions before rendering.

``` r
# Compact matrix: expanded automatically (message printed)
vismat(mat_compact,labelcex=0.5)
#> Expanding compact matrix (27 -> 28 individuals) for visualization.
```

![](relationship-matrix_files/figure-html/heatmap_compact-1.png)

#### Preserve Pedigree Order

Set `reorder = FALSE` to keep the original pedigree order instead of
re-sorting by clustering.

``` r
vismat(mat_A, reorder = FALSE, labelcex = 0.5)
```

![](relationship-matrix_files/figure-html/heatmap_no_reorder-1.png)

#### Display a Subset of Individuals

Use `ids` to focus on specific individuals.

``` r
target_ids <- rownames(as.matrix(mat_A))[1:8]
vismat(mat_A, ids = target_ids,
       main = "Relationship Heatmap — First 8 Individuals")
```

![](relationship-matrix_files/figure-html/heatmap_ids-1.png)

#### Grouping by Pedigree Column

For large populations, aggregate relationships to a group-level view
using the `by` parameter. The matrix is reduced to mean coefficients
between groups.

``` r
# Mean relationship between generations
vismat(mat_A, ped = tped, by = "Gen",
       main = "Mean Relationship Between Generations")
#> Aggregating 28 individuals into 6 groups based on 'Gen'...
```

![](relationship-matrix_files/figure-html/heatmap_group-1.png)

``` r
# Mean relationship between full-sib families
# (founders without a family assignment are excluded automatically)
vismat(mat_A, ped = tped, by = "Family",
       main = "Mean Relationship Between Full-Sib Families")
#> Note: Excluding 9 founder(s) with no family assignment: J1, O, N, F, R (and 4 more)
#> Aggregating 19 individuals into 11 groups based on 'Family'...
```

![](relationship-matrix_files/figure-html/heatmap_family-1.png)

### 4.2 Inbreeding and Kinship Histograms

The “histogram” type displays the distribution of relationship
coefficients (lower triangle) or inbreeding coefficients.

``` r
# Distribution of relationship coefficients
vismat(mat_A, type = "histogram")
```

![](relationship-matrix_files/figure-html/histogram-1.png)

## 5. Performance Considerations

Calculation and visualization of large matrices can be
resource-intensive.
[`vismat()`](https://luansheng.github.io/visPedigree/reference/vismat.md)
applies the following automatic optimizations:

| Condition                    | Behavior                                                                                                                     |
|:-----------------------------|:-----------------------------------------------------------------------------------------------------------------------------|
| Compact + `by`               | Group means are computed directly from the compact matrix (no full expansion)                                                |
| Compact, no `by`, N \> 5 000 | Uses compact representative view (labels show `ID (×n)`)                                                                     |
| Compact, no `by`, N ≤ 5 000  | Matrix is automatically expanded via [`expand_pedmat()`](https://luansheng.github.io/visPedigree/reference/expand_pedmat.md) |
| N \> 2 000                   | Hierarchical clustering (reorder) is automatically skipped                                                                   |
| N \> 500                     | Individual labels are automatically hidden                                                                                   |
| N \> 100                     | Grid lines are automatically hidden                                                                                          |

When a compact `pedmat` is used with `by`,
[`vismat()`](https://luansheng.github.io/visPedigree/reference/vismat.md)
computes the group-level mean relationship matrix algebraically from the
K×K compact matrix, including a sibling off-diagonal correction. This
avoids expanding to the full N×N matrix, making family-level or
generation-level visualization feasible even for pedigrees with hundreds
of thousands of individuals.

The example below uses `big_family_size_ped` (178 431 individuals,
compact to 2 626) and displays the mean additive relationship among
**all** full-sib families in the latest generation — a computation that
would be infeasible with full expansion.

``` r
data(big_family_size_ped)

tp_big <- tidyped(big_family_size_ped)
last_gen <- max(tp_big$Gen, na.rm = TRUE)

# Compute the compact A matrix for the entire pedigree
mat_big_compact <- pedmat(tp_big, method = "A", compact = TRUE)

# Focus on all individuals in the last generation that belong to a family
ids_last_gen <- tp_big[Gen == last_gen & !is.na(Family), Ind]

# vismat() aggregates directly from the compact matrix — no expansion needed
vismat(
       mat_big_compact,
       ped = tp_big,
       ids = ids_last_gen,
       by = "Family",
       labelcex = 0.3,
       main = paste("Mean Relationship Between All Families in Generation", last_gen)
)
#> Aggregating 37009 individuals into 106 groups based on 'Family'...
```

![](relationship-matrix_files/figure-html/large_ped_tip-1.png)

This family-level view reveals the genetic structure among all 106
families comprising 37009 individuals, computed in seconds from the
compact matrix.

------------------------------------------------------------------------

**See Also:** -
[`vignette("tidy-pedigree", package = "visPedigree")`](https://luansheng.github.io/visPedigree/articles/tidy-pedigree.md) -
[`vignette("draw-pedigree", package = "visPedigree")`](https://luansheng.github.io/visPedigree/articles/draw-pedigree.md)
