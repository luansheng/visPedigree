# Ensure the ordering columns needed for export are present

`tidyped(..., addnum = FALSE)` omits `IndNum`, `SireNum`, and `DamNum`,
and `tidyped(..., addgen = FALSE)` omits `Gen`. These columns are
optional upstream, but the export formats rely on them for
parent-before-offspring ordering. When absent they are reconstructed
from `Ind`/`Sire`/`Dam` on a copy of the object, leaving the caller's
`ped` untouched.

## Usage

``` r
.ensure_export_cols(ped)
```

## Arguments

- ped:

  A complete tidyped object.

## Value

`ped` itself (when complete) or a copy with the missing columns filled
in.
