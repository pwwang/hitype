# Prepare the data for weight training

Prepare the data for weight training

## Usage

``` r
prepare_data_for_training(
  path_to_gs,
  exprs,
  level = 1,
  scaled = FALSE,
  clusters = NULL
)
```

## Arguments

- path_to_gs:

  Path to the gene set file without weights

- exprs:

  The expression matrix, or a seurat object (rows: genes, columns:
  samples/cells)

- level:

  The level of the gene sets to train weights for if you have multiple
  levels of gene sets.

- scaled:

  Whether the expression matrix is scaled

- clusters:

  A named vector of cluster ids If `exprs` is a seurat object, this is
  ignored. The cluster ids are taken from the seurat object.

## Value

A list with the gene sets, the z matrix and the clusters
