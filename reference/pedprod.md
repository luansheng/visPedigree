# Matrix-Free Pedigree Relationship Products

Computes products of the form \\A x\\, \\A X\\, \\A^{-1} x\\, or
\\A^{-1} X\\ directly from an ordered pedigree, without constructing the
full additive relationship matrix \\A\\. This is useful for large
pedigrees where explicitly storing the dense \\n \times n\\ matrix is
infeasible.

## Usage

``` r
pedprod(ped, x, method = c("A", "Ainv"))
```

## Arguments

- ped:

  A complete pedigree accepted by
  [`tidyped`](https://luansheng.github.io/visPedigree/reference/tidyped.md),
  preferably an existing `tidyped` object.

- x:

  A numeric or logical vector, or a numeric or logical matrix. An
  unnamed vector must have one value per pedigree individual in `IndNum`
  order. A named vector may contain a subset of individuals; omitted
  individuals receive value zero. For a matrix, the same rules apply to
  its rows and `rownames`.

- method:

  Character scalar. `"A"` (default) computes products with the additive
  relationship matrix. `"Ainv"` computes products with its inverse.

## Value

For vector input, a named numeric vector in pedigree order. For matrix
input, a numeric matrix whose rows are in pedigree order; input column
names are preserved.

## Details

For an ordered pedigree, the additive relationship matrix can be written
as \$\$A = T D T^\prime,\$\$ where \\T\\ is the transmission matrix and
\\D\\ contains Mendelian sampling variances. `pedprod()` applies these
factors through backward and forward pedigree traversals. For an \\n
\times p\\ right-hand side, computation requires \\O(np)\\ time and
working memory rather than materializing the \\O(n^2)\\ relationship
matrix.

Named inputs are aligned by individual ID. Every supplied name must
occur in `ped$Ind`, and duplicate, missing, or empty names are rejected.
Missing values and non-finite input values are also rejected to avoid
propagating an invalid value through all descendants or relatives.

## References

Colleau, J. J. (2002). An indirect approach to the extensive calculation
of relationship coefficients. *Genetics Selection Evolution*, 34,
409–421.
[doi:10.1051/gse:2002015](https://doi.org/10.1051/gse%3A2002015)

Colleau, J. J., Palhiere, I., Rodriguez-Ramilo, S. T., & Legarra, A.
(2017). A fast indirect method to compute functions of genomic
relationships concerning genotyped and ungenotyped individuals, for
diversity management. *Genetics Selection Evolution*, 49, 87.
[doi:10.1186/s12711-017-0363-9](https://doi.org/10.1186/s12711-017-0363-9)

## See also

[`pedmat`](https://luansheng.github.io/visPedigree/reference/pedmat.md)
for explicitly constructing relationship matrices and
[`pedrel`](https://luansheng.github.io/visPedigree/reference/pedrel.md)
for group relationship summaries.

## Examples

``` r
tped <- tidyped(small_ped)

# Equal contributions from two candidates; unspecified individuals are zero.
contributions <- c(Z1 = 0.5, Z2 = 0.5)
relationship_to_group <- pedprod(tped, contributions)

# Average additive relationship among the weighted candidates: c' A c.
sum(contributions * relationship_to_group[names(contributions)])
#> [1] 0.7910156

# Multiple contribution schemes can be evaluated in one call.
schemes <- cbind(
  Equal = c(A = 0.5, B = 0.5),
  SireA = c(A = 1.0, B = 0.0)
)
pedprod(tped, schemes)
#>      Equal   SireA
#> A  0.50000 1.00000
#> B  0.50000 0.00000
#> F  0.00000 0.00000
#> I  0.00000 0.00000
#> J1 0.00000 0.00000
#> J2 0.00000 0.00000
#> N  0.00000 0.00000
#> O  0.00000 0.00000
#> R  0.00000 0.00000
#> C  0.50000 0.50000
#> D  0.50000 0.50000
#> E  0.50000 0.50000
#> P  0.00000 0.00000
#> Q  0.00000 0.00000
#> G  0.25000 0.25000
#> H  0.25000 0.25000
#> K  0.25000 0.25000
#> L  0.25000 0.25000
#> M  0.25000 0.25000
#> S  0.00000 0.00000
#> T  0.00000 0.00000
#> U  0.12500 0.12500
#> V  0.25000 0.25000
#> W  0.12500 0.12500
#> X  0.18750 0.18750
#> Y  0.06250 0.06250
#> Z1 0.09375 0.09375
#> Z2 0.09375 0.09375
```
