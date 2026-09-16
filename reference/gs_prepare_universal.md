# Prepare a gene set in the universal marker format

The universal marker format is a long table (one row per gene per cell
type) with `cell_type` and `gene` columns, as used by biopipen
(<https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/>). The
columns are canonicalized by
[`canonicalize_marker_cols()`](https://pwwang.github.io/hitype/reference/canonicalize_marker_cols.md)
before this is called.

## Usage

``` r
gs_prepare_universal(cm, tissue_type = NULL, weight_encoding = function(x) x)
```

## Arguments

- cm:

  The canonicalized marker table

- tissue_type:

  The tissue type to filter the markers by

- weight_encoding:

  The weight encoding function

## Value

A list with `gene_sets` and `cell_names` (NULL)
