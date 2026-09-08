# Changelog

## Version 0.0.5

- ✨ Add
  [`find_markers()`](https://pwwang.github.io/hitype/reference/find_markers.md)
  to discover marker genes from your data with dependency-light `fc`
  (default), `Seurat`, or `presto` backends, directly producing a db
  `data.frame` consumable by
  [`gs_prepare()`](https://pwwang.github.io/hitype/reference/gs_prepare.md).
- ✨ Add `method` argument to
  [`train_weights()`](https://pwwang.github.io/hitype/reference/train_weights.md)
  with 7 weight-learning backends: `uniform`, `correlation`, `lr`,
  `glmnet` (default, recommended), `rf`, `xgb`, `lrp`. `epochs` and
  `batch_size` now apply to the `lrp` method only.
- ✨ Add `cv_folds` argument for cross-validated weight estimation:
  stratified folds, per-fold weights averaged for stability. Used by
  `lr`, `glmnet`, and `lrp`.
- ✨ Add `norm` and `use_sensitivity` arguments to
  [`hitype_score()`](https://pwwang.github.io/hitype/reference/hitype_score.md)
  and
  [`hitype_score_level()`](https://pwwang.github.io/hitype/reference/hitype_score_level.md):
  - `norm = "weight"` normalizes scores by `sum(abs(weights))` (fairer
    across cell types with different marker numbers);
  - `use_sensitivity = FALSE` disables marker-sensitivity weighting —
    recommended when scoring with learned weights to avoid
    double-penalizing shared markers.
- 🔧 Statistical soundness: removed input masking in
  [`prepare_data_for_training()`](https://pwwang.github.io/hitype/reference/prepare_data_for_training.md)
  (masked zeros let the model learn the mask pattern instead of
  expression).
- 🔧 Seurat v5 compatibility: use
  `Seurat::GetAssayData(x, layer = "data")` instead of direct
  `@assays$RNA@data` / `@assays$RNA@scale.data` access.
- 🧩 New suggested dependencies: `glmnet`, `ranger`, `xgboost` (required
  per method, optional overall).

## Version 0.0.4

- ✨ Add weight_encoding argument to gs_prepare

## Version 0.0.3

- Move keras and innsight as suggested

## Version 0.0.2

- Add logo
- Reduce keras requirement to 2.11

## Version 0.0.1

- Add `make_unique` to `RunHitype`

## Version 0.0.0

- Initial release
