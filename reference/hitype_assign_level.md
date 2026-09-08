# Assign cell types based on ScType scores

Assign cell types based on ScType scores

## Usage

``` r
hitype_assign_level(
  clusters,
  scores,
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

  A matrix of cell type scores (cell_types x cells)

- threshold:

  Confidence threshold as top1/top2 score ratio. When `NULL` (default),
  no filtering is done. When a number, the rank-1 cell type of a cluster
  is marked as `<UNKNOWN>` when its score is less than `threshold` times
  the second-best score.

- top:

  The number of top cell types to assign for each cluster

- mode:

  `"cluster"` (default) aggregates scores per cluster then assigns.
  `"cell"` assigns each cell individually then reports majority vote per
  cluster.

## Value

A data from of top cell type assignments with columns: Cluster,
CellType, Score, Margin
