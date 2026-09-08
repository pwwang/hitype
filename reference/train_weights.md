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
  range = c(1, 5),
  data_split = c(0.7, 0.2, 0.1),
  epochs = 20,
  batch_size = 32,
  run_weights_on_test = TRUE,
  cv_folds = 1,
  method = c("glmnet", "lr", "rf", "xgb", "lrp", "correlation", "uniform")
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

- range:

  The range of the weights

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

## Value

A data frame with the weights, that can be used directly by
[`gs_prepare()`](https://pwwang.github.io/hitype/reference/gs_prepare.md).
