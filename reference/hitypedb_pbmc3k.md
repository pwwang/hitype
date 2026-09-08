# Gene sets with weights trained from the PBMC 3k dataset

Marker gene sets for the 9 cell types of the Seurat PBMC 3k tutorial
(`Naive CD4+ T`, `CD14+ Mono`, `Memory CD4+`, `B`, `CD8+ T`,
`FCGR3A+ Mono`, `NK`, `DC` and `Platelet`), with weights trained on the
dataset itself. It is intended for demonstration: pass it to
[`gs_prepare`](https://pwwang.github.io/hitype/reference/gs_prepare.md)
like any other marker data frame (e.g. `gs_prepare(hitypedb_pbmc3k)`).

## Usage

``` r
hitypedb_pbmc3k
```

## Format

A data frame in the marker database format used by
[`gs_prepare`](https://pwwang.github.io/hitype/reference/gs_prepare.md):
`cellName`, `geneSymbolmore1`, `geneSymbolmore2` and `level`, where
`geneSymbolmore1` carries `+`/`-` suffixes encoding the trained weights.

## Source

Trained with
[`train_weights`](https://pwwang.github.io/hitype/reference/train_weights.md)
on the PBMC 3k dataset (Seurat guided-clustering tutorial, resolution
0.5).
