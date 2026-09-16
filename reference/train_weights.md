# Train weights for the markers

Train weights for the markers

## Usage

``` r
train_weights(
  path_to_gs,
  exprs,
  level = 1,
  scaled = FALSE,
  clusters = NULL,
  range = c(-5, -1, 1, 5),
  data_split = c(0.7, 0.2, 0.1),
  epochs = 20,
  batch_size = 32,
  run_weights_on_test = TRUE,
  cv_folds = 1,
  method = c("glmnet", "lr", "rf", "xgb", "lrp", "correlation", "uniform"),
  format = c("universal", "db"),
  drop_zero = TRUE,
  pos_only = FALSE,
  return_models = FALSE,
  seed = 8525
)
```

## Arguments

- path_to_gs:

  Path to the gene set file without weights. The training cell types are
  the cell types of the cells (see `clusters`): a cell type in the
  marker file with cells of that type in the data is trained on its own
  markers in the file. A cell type in the marker file with no cells of
  that type in the data is ignored with a warning, and its markers are
  pooled for the cell types in the data that are not covered by the
  marker file (each of them is trained on the pooled markers). Cell
  types in the data that are neither covered by the marker file nor by
  pooled markers are not trained (with a warning).

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

- range:

  The quantization range for the db output only (`format = "db"`): the
  4-element form `c(-low, -high, low, high)` (default `c(-5, -1, 1, 5)`)
  scales the positive weights into `[low, high]` and the negative ones
  into `[-high, -low]` before rounding them to integers, so genuinely
  positive markers always stay positive. The universal output (default)
  keeps the raw trained weights as-is.

- data_split:

  A vector of fractions for training, validation and testing. If only
  two fractions are provided, no testing set will be used.

- epochs:

  The number of epochs to train (lrp method only)

- batch_size:

  The batch size (lrp method only)

- run_weights_on_test:

  Whether to run the weights on the test set. Requires that `data_split`
  has three elements.

- cv_folds:

  Number of cross-validation folds for weight estimation. When \> 1,
  weights are averaged across folds for stability. Default is 1 (no
  cross-validation). Used by lr, glmnet, and lrp methods.

- method:

  The weight learning method. One of:

  "uniform"

  :   All markers get equal weight (= 1). Fast baseline.

  "correlation"

  :   Pearson correlation between each marker and the binary cluster
      indicator. Simple, interpretable.

  "lr"

  :   Logistic regression coefficients (one-vs-rest). Classic ML
      approach. Use cv_folds for stability.

  "glmnet"

  :   Sparse logistic regression with elastic net penalty (alpha = 0.5).
      Automatically zeros out uninformative markers. Recommended for
      most users. Requires the glmnet package.

  "rf"

  :   Random forest permutation importance. Captures non-linear marker
      interactions. Requires the ranger package.

  "xgb"

  :   XGBoost gain-based feature importance. State-of-the-art tree
      method. Requires the xgboost package.

  "lrp"

  :   Neural network + Layer-wise Relevance Propagation. Deep learning
      approach. Requires the keras and innsight packages.

- format:

  The format of the output data frame. One of `"universal"` (default) or
  `"db"` (the hitype/ScType wide format).

- drop_zero:

  Whether to drop markers whose trained weight is exactly zero from the
  output (default `TRUE`). Such markers carry no learned direction —
  with the glmnet method (coefficients at `lambda.1se`) most candidate
  markers are zeroed. Dropping them makes zero mean "no direction": the
  marker is absent from the cell type's list. Pass `FALSE` to keep every
  candidate marker (with `format = "db"` they are then encoded as `*`).

- pos_only:

  Whether to keep only the markers with a positive trained weight in the
  output (default `FALSE`). A marker can be trained with a negative
  weight (anti-correlating with the cell type, e.g. by the
  `"correlation"` method); such rows are reported with
  `direction = "negative"` and are dropped when `pos_only = TRUE`, so
  the returned table contains markers overexpressed in each cell type
  only.

- return_models:

  Whether to also return the fitted per-cell-type prediction models
  (default `FALSE`). With `FALSE` (default) the return value is
  unchanged: a data frame with the weights. With `TRUE`, a list is
  returned instead:

  `weights`

  :   The usual weights data frame (universal or db format).

  `models`

  :   A model bundle list that
      [`hitype_score_models()`](https://pwwang.github.io/hitype/reference/hitype_score_models.md)
      can score new cells with:

      `method`

      :   The weight learning method.

      `level`

      :   The level of the gene sets.

      `features`

      :   The ordered gene list the models were fit on.

      `center`, `scale`

      :   Numeric vectors named by `features` with the per-gene
          centering/scaling that was applied to the training matrix, so
          scoring can reproduce it exactly.

      `coefs`

      :   A named list, one entry per trained cell type, each a named
          numeric coefficient vector of the linear predictor including
          the intercept as `"(Intercept)"`. `NULL` for methods without
          linear models.

  `labels`

  :   The trained cell types (names of `coefs`).

  Model bundles are only produced for the linear-model methods
  `"glmnet"` and `"lr"`; with any other method `coefs` is `NULL`,
  `labels` is `NULL`, and
  [`hitype_score_models()`](https://pwwang.github.io/hitype/reference/hitype_score_models.md)
  errors.

- seed:

  Random seed for reproducibility

## Value

By default, a data frame with the weights in the universal marker format
(default) or the db format (`format = "db"`), that can be used directly
by
[`gs_prepare()`](https://pwwang.github.io/hitype/reference/gs_prepare.md).
With `return_models = TRUE`, a list with the weights data frame, the
model bundle and the trained cell type labels (see `return_models`).
