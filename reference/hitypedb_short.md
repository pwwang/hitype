# Gene sets from ScType (short version)

The short marker database from
[ScType](https://github.com/IanevskiAleksandr/sc-type), bundled for
convenience. It can be passed directly to
[`gs_prepare`](https://pwwang.github.io/hitype/reference/gs_prepare.md):

## Usage

``` r
hitypedb_short
```

## Format

A data frame with columns commonly used by
[`gs_prepare`](https://pwwang.github.io/hitype/reference/gs_prepare.md):
`cellName`, `tissueType`, `geneSymbolmore1`, `geneSymbolmore2` and
`level`.

## Source

<https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx>

## Details

`gs <- gs_prepare(hitypedb_short, tissue_type = "Immune system")`
