# Run compiled weights on test data

Run compiled weights on test data

## Usage

``` r
run_weights_on_test_data(db, exprs, clusters, scaled, test_data_idx)
```

## Arguments

- db:

  The weights data frame

- exprs:

  The expression matrix, or a seurat object (rows: genes, columns:
  samples/cells)

- clusters:

  A named vector of cluster ids If `exprs` is a seurat object, this is
  ignored. The cluster ids are taken from the seurat object.

- scaled:

  Whether the expression matrix is scaled

- test_data_idx:

  The row names of the test data
