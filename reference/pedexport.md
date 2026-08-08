# Export a pedigree to a breeding software format

Converts a `tidyped` pedigree into the tabular format expected by common
animal/plant breeding software packages and either returns the result as
a `data.table` or writes it to a plain-text file.

## Usage

``` r
pedexport(
  ped,
  software = c("blupf90", "asreml", "echidna", "wombat", "mtdfreml", "dmu", "numeric",
    "sommer"),
  file = NULL,
  sep = " ",
  header = NULL,
  missing = NULL
)
```

## Arguments

- ped:

  A complete `tidyped` object (see
  [`tidyped`](https://luansheng.github.io/visPedigree/reference/tidyped.md)).
  The pedigree must be structurally complete: every non-missing `Sire`
  and `Dam` identifier must be present in `Ind`.

- software:

  Character scalar specifying the target format. One of:

  `"blupf90"`

  :   Three integer columns (`IndNum`, `SireNum`, `DamNum`), no header,
      missing parents encoded as `0`. Rows are sorted by `IndNum` so
      parent rows always precede offspring rows. Compatible with
      BLUPF90, REMLF90, GIBBSF90, and related programs.

  `"asreml"`

  :   Three character columns (`animal`, `sire`, `dam`), with a one-line
      header by default, missing parents encoded as `"0"`. Rows are
      sorted by generation (`Gen`) so founders appear first. Compatible
      with ASReml-R and ASReml-SA (character IDs require the `!ALPHA`
      qualifier; the header line requires the `!SKIP 1` qualifier).

  `"echidna"`

  :   Identical character layout and defaults to `"asreml"`: columns
      `animal`, `sire`, and `dam`, a header by default, missing parents
      encoded as `"0"`, and rows sorted by generation.

  `"wombat"`

  :   Three character columns (`animal`, `sire`, `dam`), no header,
      missing parents encoded as `"0"`. Wombat accepts alphanumeric IDs
      and recodes them internally. Rows are sorted so parents precede
      offspring.

  `"mtdfreml"`

  :   Identical integer layout to `"blupf90"`. Compatible with MTDFREML
      and MTGSAM.

  `"dmu"`

  :   Identical integer layout to `"blupf90"`. Compatible with DMU.

  `"numeric"`

  :   Integer layout with an optional header (`IndNum SireNum DamNum`).
      A portable, software-agnostic choice when the target program is
      not listed above.

  `"sommer"`

  :   Three character columns (`ID`, `Sire`, `Dam`), missing parents
      coded as `NA`. Rows are sorted so parents precede offspring.
      Returned as a data.table ready to pair with a base R relationship
      matrix from
      [`pedmat`](https://luansheng.github.io/visPedigree/reference/pedmat.md)`(..., sparse = FALSE)`
      for use in `sommer::mmer(..., random = ~vsr(animal, Gu = A))`.

- file:

  Character scalar. Path to the output file. When `NULL` (default) the
  formatted `data.table` is returned invisibly and no file is written.
  File output is not supported for `"sommer"`, which is intended for
  direct use as an R `data.table`.

- sep:

  Character scalar. Field separator used when writing to `file`.
  Defaults to a single space (`" "`), which every file-based format
  accepts. BLUPF90 requires a space; WOMBAT, MTDFREML, and DMU accept
  spaces or TABs; ASReml and Echidna accept spaces, TABs, or commas; and
  the generic `"numeric"` format accepts any single-byte, non-empty
  separator supported by
  [`fwrite`](https://rdrr.io/pkg/data.table/man/fwrite.html).
  ASReml/Echidna comma-delimited output requires a `.csv` file extension
  or the `!CSV` qualifier.

- header:

  Logical scalar or `NULL`. Whether to include a column header line.
  `NULL` (default) uses the software-specific default: `TRUE` for
  `"asreml"`, `"echidna"`, and `"numeric"`, `FALSE` for `"blupf90"`,
  `"wombat"`, `"mtdfreml"` and `"dmu"`. Ignored for `"sommer"`. Note for
  ASReml and Echidna: a header line is read as data unless the `!SKIP 1`
  qualifier is used in the command file.

- missing:

  Character or integer scalar. Symbol for missing parents. `NULL`
  (default) uses the software-specific default: `0L` for numeric
  formats, `"0"` for `"asreml"`, `"echidna"`, and `"wombat"`. Numeric
  formats (`"blupf90"`, `"mtdfreml"`, `"dmu"`, `"numeric"`) require a
  single integer value; `"asreml"`, `"echidna"`, and `"wombat"` accept a
  character value (numeric values are converted to character). Ignored
  for `"sommer"`, which always codes missing parents as `NA`.

## Value

A `data.table` in the target format, returned invisibly. Numeric formats
carry an `xref` attribute mapping each numeric ID back to its original
character ID (columns `IndNum` and `Ind`), and when `file` is given the
mapping is also written to `paste0(file, ".xref")` without a header,
mirroring RENUMF90's `_XrefID` file. For file-based formats, when `file`
is not `NULL`, the table is also written to that path and a message
reports the number of individuals written.

## Details

**Column semantics by format:**

|  |  |  |  |
|----|----|----|----|
| **software** | **Col 1** | **Col 2** | **Col 3** |
| blupf90 / mtdfreml / dmu / numeric | `IndNum` (integer) | `SireNum` (integer) | `DamNum` (integer) |
| asreml / echidna / wombat | `animal` (character) | `sire` (character) | `dam` (character) |
| sommer | `ID` (character) | `Sire` (character) | `Dam` (character) |

**Software format requirements:**

All file-based formats use unquoted, free-format plain text. The table
below summarises what each program expects, and the defaults
`pedexport()` uses to meet those expectations.

|  |  |  |  |  |  |
|----|----|----|----|----|----|
| **software** | **header** | **separator** | **missing** | **IDs** | **notes** |
| blupf90 | no | spaces only | `0` | integer | TAB separators rejected |
| wombat | no | space / TAB | `"0"` | character | alphanumeric accepted, recoded internally |
| mtdfreml | no | space / TAB | `0` | integer |  |
| dmu | no | space / TAB | `0` | integer |  |
| asreml | yes | space / TAB / comma | `"0"` | char or integer | character IDs need `!ALPHA`; header needs `!SKIP 1`; comma needs `.csv` or `!CSV` |
| echidna | yes | space / TAB / comma | `"0"` | char or integer | same pedigree layout and defaults as ASReml |
| sommer | n/a | n/a | `NA` | character | returned as a data.table, not a file |
| numeric | optional | any | `0` | integer | generic software-agnostic layout |

**Row order:**

All numeric formats sort by `IndNum` ascending, which guarantees that
parent rows appear before offspring rows (parents always receive a
smaller integer index after
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
topological sorting). The `"asreml"`, `"echidna"`, `"wombat"` and
`"sommer"` formats sort by `Gen` ascending for the same reason.

**Optional
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
columns:**

`tidyped(..., addnum = FALSE)` omits `IndNum`/`SireNum`/ `DamNum` and
`tidyped(..., addgen = FALSE)` omits `Gen`. When these columns are
absent, `pedexport()` reconstructs them from `Ind`/`Sire`/`Dam`. Because
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
always topologically sorts rows, parents still precede offspring in the
output.

**Identifiers in written files:**

Output is deliberately unquoted because the supported breeding programs
read pedigree files as free-format text. Therefore, identifiers and
character missing-parent symbols must not contain whitespace. They also
must not contain the selected non-whitespace separator. Numeric formats
apply the same restriction to original identifiers written to the
companion `.xref` file.

**Sommer compatibility:**

The `Gu` argument of `sommer::vsr()` requires a base R matrix. Because
[`pedmat()`](https://luansheng.github.io/visPedigree/reference/pedmat.md)
defaults to `sparse = TRUE` and returns a Matrix object for
`method = "A"`, use `pedmat(ped, method = "A", sparse = FALSE)` or
explicitly convert the result with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html). The matrix row and
column names match the `ID` column returned by
`pedexport(..., software = "sommer")`.

**Relationship to
[`tidyped()`](https://luansheng.github.io/visPedigree/reference/tidyped.md):**

`pedexport()` is a post-processing step that re-encodes an already
validated pedigree. Pass the output of
[`tidyped`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
directly to avoid redundant validation.

## See also

[`tidyped`](https://luansheng.github.io/visPedigree/reference/tidyped.md)
to prepare the pedigree,
[`pedmat`](https://luansheng.github.io/visPedigree/reference/pedmat.md)
to compute relationship matrices.

## Examples

``` r
library(visPedigree)

tp <- tidyped(small_ped)

# Return a data.table (no file written)
out_blupf90 <- pedexport(tp, software = "blupf90")
head(out_blupf90)
#>    IndNum SireNum DamNum
#>     <int>   <int>  <int>
#> 1:      1       0      0
#> 2:      2       0      0
#> 3:      3       0      0
#> 4:      4       0      0
#> 5:      5       0      0
#> 6:      6       0      0

# ASReml format with character IDs and header
out_asreml <- pedexport(tp, software = "asreml")
head(out_asreml)
#>    animal   sire    dam
#>    <char> <char> <char>
#> 1:      A      0      0
#> 2:      B      0      0
#> 3:      F      0      0
#> 4:      I      0      0
#> 5:     J1      0      0
#> 6:     J2      0      0

# Echidna uses the same pedigree layout as ASReml
out_echidna <- pedexport(tp, software = "echidna")
identical(out_echidna, out_asreml)
#> [1] TRUE

# Numeric formats carry the ID mapping back to character IDs
out_dmu <- pedexport(tp, software = "dmu")
head(attr(out_dmu, "xref"))
#>    IndNum    Ind
#>     <int> <char>
#> 1:      1      A
#> 2:      2      B
#> 3:      3      F
#> 4:      4      I
#> 5:      5     J1
#> 6:      6     J2

# sommer: character pedigree + visPedigree's own relationship matrix.
# pedmat() rownames match pedexport(..., software = "sommer")$ID, so the
# dense matrix is ready for sommer::mmer(..., Gu = A).
out_sommer <- pedexport(tp, software = "sommer")
head(out_sommer)
#>        ID   Sire    Dam
#>    <char> <char> <char>
#> 1:      A   <NA>   <NA>
#> 2:      B   <NA>   <NA>
#> 3:      F   <NA>   <NA>
#> 4:      I   <NA>   <NA>
#> 5:     J1   <NA>   <NA>
#> 6:     J2   <NA>   <NA>
A <- pedmat(tp, method = "A", sparse = FALSE)
identical(rownames(A), out_sommer$ID)
#> [1] TRUE

# \donttest{
# Write to a file (numeric formats also write <file>.xref)
tmp <- tempfile(fileext = ".txt")
pedexport(tp, software = "blupf90", file = tmp)
#> Written ID mapping to: /tmp/RtmpgWc7jQ/file1c4551c85467.txt.xref
#> Written 28 individuals to: /tmp/RtmpgWc7jQ/file1c4551c85467.txt
readLines(tmp, n = 5)
#> [1] "1 0 0" "2 0 0" "3 0 0" "4 0 0" "5 0 0"
readLines(paste0(tmp, ".xref"), n = 5)
#> [1] "1 A"  "2 B"  "3 F"  "4 I"  "5 J1"
# }
```
