# Prepare gene sets

## Using gene sets from [ScType](https://github.com/IanevskiAleksandr/sc-type)

`hitype` is designed to be compatible with
[ScType](https://github.com/IanevskiAleksandr/sc-type). So the gene sets
provided by [ScType](https://github.com/IanevskiAleksandr/sc-type) can
be used directly.

``` r
library(hitype)

gs <- gs_prepare(
  "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx"
)
gs$gene_sets[[1]]$`Acinar cells`
#> $markers
#>  [1] "CTRB1"    "KLK1"     "RBPJL"    "PTF1A"    "CELA3A"   "PRSS1"   
#>  [7] "SPINK1"   "ZG16"     "CEL"      "CELA2A"   "CPB1"     "CELA1"   
#> [13] "RNASE1"   "AMY2B"    "CPA2"     "CPA1"     "CELA3B"   "PNLIP"   
#> [19] "CTRB2"    "PLA2G1B"  "PRSS2"    "CLPS"     "REG1A"    "SYCN"    
#> [25] "PNLIPRP1" "CTRC"     "REG3A"    "PRSS3"    "REG1B"    "CFB"     
#> [31] "GDF15"    "MUC1"     "C15orf48" "AKR1C3"   "OLFM4"    "GSTA1"   
#> [37] "LGALS2"   "PDZK1IP1" "RARRES2"  "CXCL17"   "GSTA2"    "ANPEP"   
#> [43] "LYZ"      "ANGPTL4"  "ALDOB"   
#> 
#> $weights
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [39] 1 1 1 1 1 1 1

# gs <- gs_prepare(
#  "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_long.xlsx"
# )
```

## Preparing your own gene sets

You can also prepare your own gene sets. The gene sets should be a
`data.frame` or a tab-delimited file with the following columns:

- `tissueType` (optional): Tissue type. One can pass `tissue_type` to
  [`gs_prepare()`](https://pwwang.github.io/hitype/reference/gs_prepare.md)
  to filter gene sets by tissue type.
- `cellName` (required): Cell type name.
- `geneSymbolmore1` (required): Markers for the cell type. Multiple
  markers should be separated by `,`. One can use a suffix `+` to
  indicate that the marker is a positive marker, and `-` to indicate
  that the marker is a negative marker. Multiple `+`’s are allowed to
  indicate that the marker is a strong positive marker. Multiple `-`’s
  are allowed to indicate that the marker is a strong negative marker.
  For example, `CD3E+++` indicates that `CD3E` is a strong positive
  marker, and `CD14---` indicates that `CD14` is a strong negative
  marker.
- `geneSymbolmore2` (required): Negative markers for the cell type. This
  column can be empty strings. This is kept for compatibility with
  [ScType](https://github.com/IanevskiAleksandr/sc-type). No suffixes
  are allowed in this column. A marker in this column is the same as a
  marker in the previous column with a suffix `-`.
- `level` (optional): Cell type level. The level of the `cellName`. It
  must start from 1 and increase by 1. The final cell type consists of
  one `cellName` from each level.
- `nextLevels`: Indication of possible `cellName`s at the next level.
  Levels should be separated by `;`. The `cellName` at each level should
  be separated by `,`.
  - The format is `cellName1,cellName2;cellName3,cellName4;...` for each
    `cellName` at the current level.
  - For example, for `CD4` at level 1, `nextLevels` `Naive;Activated`
    indicates that `CD4` should be followed by `Naive` at level 2 and
    `Activated` at level 3.
  - If the `nextLevels` is also specified for `Naive` at level 2, then
    the `cellName` at level 3 for `CD4` will be the intersection of
    `Activated` and the `cellName`s at level 3 for `Naive` limited by
    `nextLevels` for `Naive`.
  - If a level in `nextLevels` is empty, then the `cellName`s at that
    level will be all `cellName`s at the next level limited by
    `nextLevels` of that level.
  - `!` can be used to indicate that the `cellName`s at the next level
    should be excluded. For example, `!Naive` indicates that `Naive`
    should be excluded from the `cellName`s at the next level.
  - If `!` is the only character in a level, then all `cellName`s at
    that level will be excluded.

Any `data.frame` of markers (e.g. the output of
[`train_weights()`](https://pwwang.github.io/hitype/reference/train_weights.md))
can be passed directly to
[`gs_prepare()`](https://pwwang.github.io/hitype/reference/gs_prepare.md):

``` r
markers <- data.frame(
    cellName = c("CD4 T", "CD8 T"),
    geneSymbolmore1 = c("CD3E,CD4,IL7R", "CD3E,CD8A,CD8B"),
    geneSymbolmore2 = c("", "")
)
gs <- gs_prepare(markers)
```
