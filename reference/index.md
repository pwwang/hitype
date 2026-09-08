# Package index

## Cell type scoring and assignment with Hitype

- [`gs_prepare()`](https://pwwang.github.io/hitype/reference/gs_prepare.md)
  : Prepare gene sets for hitype
- [`hitype_score()`](https://pwwang.github.io/hitype/reference/hitype_score.md)
  : Calculate cell type scores
- [`hitype_assign()`](https://pwwang.github.io/hitype/reference/hitype_assign.md)
  : Generate scores for cell types for each level
- [`RunHitype()`](https://pwwang.github.io/hitype/reference/RunHitype.md)
  : Run hitype_assign for a Seurat object

## Marker weight training

- [`train_weights()`](https://pwwang.github.io/hitype/reference/train_weights.md)
  : Train weights for the markers
- [`find_markers()`](https://pwwang.github.io/hitype/reference/find_markers.md)
  : Find marker genes for cell types

## Summarizing Hitype results

- [`summary(`*`<hitype_result>`*`)`](https://pwwang.github.io/hitype/reference/summary.hitype_result.md)
  : Summarize the hitype_result object
- [`print(`*`<hitype_result>`*`)`](https://pwwang.github.io/hitype/reference/print.hitype_result.md)
  : Print the summary of the hitype_result object

## Bundled marker databases

- [`hitypedb_full`](https://pwwang.github.io/hitype/reference/hitypedb_full.md)
  : Gene sets from ScType (full version)
- [`hitypedb_short`](https://pwwang.github.io/hitype/reference/hitypedb_short.md)
  : Gene sets from ScType (short version)
- [`hitypedb_pbmc3k`](https://pwwang.github.io/hitype/reference/hitypedb_pbmc3k.md)
  : Gene sets with weights trained from the PBMC 3k dataset
