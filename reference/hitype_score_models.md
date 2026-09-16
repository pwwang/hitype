# Score cells with the models trained by train_weights()

Scores each cell with the linear predictors of the per-cell-type models
persisted by
[`train_weights()`](https://pwwang.github.io/hitype/reference/train_weights.md)
with `return_models = TRUE` (methods `"glmnet"` and `"lr"` only): for a
cell with expression `x` of the model features, the score of cell type
`t` is
`eta_t = (Intercept)_t + sum_f coef_tf * (x_f - center_f) / scale_f`,
where `center`/`scale` are the per-gene centering/scaling recorded when
the models were trained. Genes of the model features missing from
`exprs` contribute 0; extra genes of `exprs` are ignored.

## Usage

``` r
hitype_score_models(exprs, models, margin = 0)
```

## Arguments

- exprs:

  Input scRNA-seq expression matrix (genes x cells, the same convention
  as
  [`hitype_score()`](https://pwwang.github.io/hitype/reference/hitype_score.md)).

- models:

  The model bundle returned by
  [`train_weights()`](https://pwwang.github.io/hitype/reference/train_weights.md)
  with `return_models = TRUE` (the `models` element).

- margin:

  Cells whose top-minus-second score (`margins`) is below `margin` are
  assigned `"Unknown"` instead of their top cell type. `0` (default)
  assigns every cell.

## Value

A list with:

- `scores`:

  A matrix (cells x cell types) of the linear predictors, with the cell
  types in the order of the model bundle.

- `assignments`:

  A named vector with the top-scoring cell type of every cell (or
  `"Unknown"` for cells below the `margin`).

- `margins`:

  A named numeric vector with the top-minus-second score of every cell.

## Details

The centering/scaling is folded into the coefficients, so a sparse
`exprs` matrix is never densified by the `(x - center) / scale` shift:
`eta_t = (Intercept)_t - sum_f b_tf * center_f + sum_f b_tf * x_f` with
`b_tf = coef_tf / scale_f`.
