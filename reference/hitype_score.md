# Calculate cell type scores

Calculate cell type scores

## Usage

``` r
hitype_score(exprs, gs, scaled = FALSE, norm = "sqrt", use_sensitivity = TRUE)
```

## Arguments

- exprs:

  Input scRNA-seq data matrix (rownames - gene names, colnames - cell
  names)

- gs:

  The gene sets prepared by
  [`gs_prepare`](https://pwwang.github.io/hitype/reference/gs_prepare.md)
  The `gene_sets` is used. One could also pass `gs$gene_sets` directly.

- scaled:

  Whether the input data is scaled or not

- norm:

  Normalization method for scoring. `"sqrt"` (default) divides by
  `sqrt(n)`; `"weight"` divides by `sum(abs(weights))`.

- use_sensitivity:

  Whether to multiply expression z-scores by marker sensitivity scores.
  Set to `FALSE` when using learned weights to avoid double-penalizing
  shared markers. Default is `TRUE`.

## Value

A list of matrices of cell type scores for each level (rownames - cell
types, colnames - cell names)

## Author

Matt Mulvahill, Panwen Wang
