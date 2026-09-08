# Gene sets from ScType (full version)

The full marker database from
[ScType](https://github.com/IanevskiAleksandr/sc-type), bundled for
convenience. It can be passed directly to
[`gs_prepare`](https://pwwang.github.io/hitype/reference/gs_prepare.md),
e.g. `gs_prepare(hitypedb_full, tissue_type = "Immune system")`.

## Usage

``` r
hitypedb_full
```

## Format

A data frame with columns commonly used by
[`gs_prepare`](https://pwwang.github.io/hitype/reference/gs_prepare.md):
`cellName`, `tissueType`, `geneSymbolmore1`, `geneSymbolmore2` and
`level`.

## Source

<https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx>
