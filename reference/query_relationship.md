# Query Relationship Coefficients from a Pedigree Matrix

Retrieves relationship coefficients between individuals from a pedmat
object. For compact matrices, automatically handles lookup of merged
full-siblings.

## Usage

``` r
query_relationship(x, id1, id2 = NULL)
```

## Arguments

- x:

  A pedmat object created by
  [`pedmat`](https://luansheng.github.io/visPedigree/reference/pedmat.md).

- id1:

  Character, first individual ID.

- id2:

  Character, second individual ID. If `NULL`, returns the entire row of
  relationships for `id1`.

## Value

- If `id2` is provided: numeric value (relationship coefficient)

- If `id2` is `NULL`: named numeric vector (id1's row)

- Returns `NA` if individual not found

## Details

For compact matrices (`compact = TRUE`), this function automatically
maps individuals to their family representatives. For methods A, D, and
AA, it can compute the correct relationship even between merged
full-siblings using the formula:

- **A**: \\a\_{ij} = 0.5 \* (a\_{is} + a\_{id})\\ where s, d are parents

- **D**: \\d\_{ij} = a\_{ij}^2\\ (for full-sibs in same family)

- **AA**: \\aa\_{ij} = a\_{ij}^2\\

## Note

Inverse matrices (Ainv, Dinv, AAinv) are **not supported** because
inverse matrix elements do not represent meaningful relationship
coefficients.

## See also

[`pedmat`](https://luansheng.github.io/visPedigree/reference/pedmat.md),
[`expand_pedmat`](https://luansheng.github.io/visPedigree/reference/expand_pedmat.md)

## Examples

``` r
tped <- tidyped(small_ped)
A <- pedmat(tped, method = "A", compact = TRUE)

# Query specific pair
query_relationship(A, "A", "B")
#> [1] 0

# Query merged full-siblings (works with compact)
query_relationship(A, "Z1", "Z2")
#> [1] 0.5507812

# Get all relationships for one individual
query_relationship(A, "A")
#>       A       B       F       I      J1      J2       N       O       R       C 
#> 1.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.50000 
#>       D       E       P       Q       G       H       K       L       M       S 
#> 0.50000 0.50000 0.00000 0.00000 0.25000 0.25000 0.25000 0.25000 0.25000 0.00000 
#>       T       U       V       W       X       Y      Z1 
#> 0.00000 0.12500 0.25000 0.12500 0.18750 0.06250 0.09375 
```
