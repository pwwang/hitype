# Generate scores for cell types for each level

Generate scores for cell types for each level

## Usage

``` r
hitype_assign(
  clusters,
  scores,
  gs = NULL,
  fallback = "Unknown",
  threshold = NULL,
  top = 10,
  mode = c("cluster", "cell")
)
```

## Arguments

- clusters:

  A named vector of original cluster assignments (names - cell names,
  values - cluster assignments)

- scores:

  A list of matrices of cell type scores for each level

- gs:

  The gene sets prepared by
  [`gs_prepare`](https://pwwang.github.io/hitype/reference/gs_prepare.md)
  The `cell_names` is actually used. One could also pass `gs$cell_names`
  directly.

- fallback:

  A fallback cell type if no cell type is assigned

- threshold:

  Confidence threshold as top1/top2 score ratio. `NULL` (default) means
  no confidence filtering. A number marks the rank-1 cell type of a
  cluster as `<UNKNOWN>` when its score is less than `threshold` times
  the second-best score.

- top:

  The number of top cell types to assign for each cluster in the result.

- mode:

  `"cluster"` (default) aggregates scores per cluster then assigns the
  top cell type. `"cell"` assigns each cell individually then reports
  majority vote per cluster.

## Value

A dataframe with columns: `Level`, `Cluster`, `CellType`, `Score`, and
`Margin`. For each level and cluster, the top cell types are returned.
You can use
[`summary.hitype_result`](https://pwwang.github.io/hitype/reference/summary.hitype_result.md)
to print the combination of cell types for each cluster.
