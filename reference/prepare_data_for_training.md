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

  The cell types of the cells. When `exprs` is a Seurat object, it can
  be `NULL` (default) to take the cell types from
  [`Seurat::Idents()`](https://satijalab.github.io/seurat-object/reference/Idents.html),
  or a column name in the `meta.data` of the Seurat object that holds
  the cell type of each cell. When `exprs` is a matrix, it should be a
  named vector of cell types (names - cell names).

## Value

A list with the gene sets, the z matrix and the clusters
