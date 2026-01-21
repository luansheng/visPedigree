# Compact pedigree by merging full siblings for matrix calculation

This internal function identifies full siblings (individuals sharing the
same sire and dam) and selects one representative per family. This can
dramatically reduce memory requirements when calculating relationship
matrices for pedigrees with large full-sibling families.

## Usage

``` r
compact_ped_for_matrix(ped)
```

## Arguments

- ped:

  A tidyped object or pedigree data.

## Value

A list containing:

- `ped_compact`: The compacted pedigree (tidyped object)

- `compact_map`: Individual-level mapping table

- `family_summary`: Family-level summary table

- `compact_stats`: Compression statistics
