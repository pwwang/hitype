# Uniform weights (baseline)

Every marker of a cell type gets the same weight, so the cell type is
defined by its own markers only.

## Usage

``` r
train_uniform(data, clusters, class_markers = NULL)
```

## Arguments

- data:

  The prepared training data, see
  [`prepare_data_for_training()`](https://pwwang.github.io/hitype/reference/prepare_data_for_training.md)

- clusters:

  The cell type of each cell of `data$z`

- class_markers:

  A named list of the marker genes of each cell type (names are the cell
  types, values the gene vectors). Each cell type then reports its own
  markers instead of every candidate feature. `NULL` (default) keeps the
  previous behaviour of reporting every candidate feature for every cell
  type. See
  [`class_features()`](https://pwwang.github.io/hitype/reference/class_features.md).
