# Calculate ScType scores and assign cell types for one level

Calculate ScType scores and assign cell types for one level

## Usage

``` r
hitype_score_level(z, gs_level, norm = "sqrt", use_sensitivity = TRUE)
```

## Arguments

- z:

  Z-scaled expression matrix (rownames - gene names, colnames - cell
  names)

- gs_level:

  One level of gene sets prepared by
  [`gs_prepare`](https://pwwang.github.io/hitype/reference/gs_prepare.md)

- norm:

  Normalization method for scoring. `"sqrt"` (default) divides by
  `sqrt(n)`; `"weight"` divides by `sum(abs(weights))`.

- use_sensitivity:

  Whether to multiply expression z-scores by marker sensitivity scores.
  Set to `FALSE` when using learned weights to avoid double-penalizing
  shared markers. Default is `TRUE`.

## Value

A matrix of cell type scores for each cell (rownames - cell types,
colnames - cell names)
