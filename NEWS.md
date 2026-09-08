## Version 0.0.7

- ✨ `train_weights()` now derives the training cell types from the data
  instead of the marker file. With a Seurat `exprs`, `clusters` can be
  `NULL` (default, using `Seurat::Idents()`) or a column name in the
  `meta.data`; with a matrix, it stays a named per-cell vector. A cell
  type in the marker file with no cells of that type in the data is
  ignored (with a warning), and its markers are pooled for the cell types
  in the data that are not covered by the marker file, so those cell types
  can still be trained. No cells are dropped from the training data, and
  the returned weights cover only the cell types present in the data.

- 🐛 Fix `train_weights()` failing on Seurat input (`exprs` as a Seurat
  object): clusters were converted to integer codes while `data$clusters`
  stayed a factor, so the per-method `clusters == ct` comparisons
  (against `unique(data$clusters)`) were all-FALSE, every supervised fit
  errored inside a `tryCatch` and was silently swallowed, and all trained
  weights came out as the rescale midpoint (e.g. uniform `3` for the
  default `range = c(1, 5)`). Clusters now stay as character cell-type
  labels throughout, which `compile_weights()` also needs to match
  `output_node` against the gene-set cell-type names.

## Version 0.0.6

- ✨ `gs_prepare()` now auto-detects the **universal marker format** (a
  long table with `cell_type` and `gene` columns, shared with biopipen's
  `CellTypeAnnotation`) in addition to the native wide ScType-style
  format. Column aliases are matched case-insensitively
  (`celltype`/`type` → `cell_type`; `marker`/`gene_symbol` → `gene`;
  `sign` → `direction`; `tissueType` → `tissue`). Optional columns:
  `direction` (positive/negative, aliases `pos`/`neg`/`+`/`-`), `weight`
  (numeric; positive markers get `abs(weight)`, negative markers get
  `-abs(weight)`, signed weights are used as-is when `direction` is
  missing), `tissue` (filtered by `tissue_type`), and `level`. Marker
  files can be `.txt`/`.tsv`, `.csv`, `.xlsx`, `.rds`, or `.qs`/`.qs2`
  (the latter only when the `qs`/`qs2` package is installed).
- ✨ `train_weights()` and `find_markers()` gain a `format` argument
  (`"universal"` by default, `"db"` for the legacy wide format), and so
  does the internal `compile_weights()`. The default universal output
  carries exact numeric weights — the legacy db format's weight-suffix
  magnitude shift no longer applies to the default flow, including
  `run_weights_on_test_data`.

## Version 0.0.5

- ✨ Add `find_markers()` to discover marker genes from your data with
  dependency-light `fc` (default), `Seurat`, or `presto` backends, directly
  producing a db `data.frame` consumable by `gs_prepare()`.
- ✨ Add `method` argument to `train_weights()` with 7 weight-learning backends:
  `uniform`, `correlation`, `lr`, `glmnet` (default, recommended), `rf`, `xgb`, `lrp`.
  `epochs` and `batch_size` now apply to the `lrp` method only.
- ✨ Add `cv_folds` argument for cross-validated weight estimation: stratified
  folds, per-fold weights averaged for stability. Used by `lr`, `glmnet`, and `lrp`.
- ✨ Add `norm` and `use_sensitivity` arguments to `hitype_score()` and
  `hitype_score_level()`:
  - `norm = "weight"` normalizes scores by `sum(abs(weights))` (fairer across
    cell types with different marker numbers);
  - `use_sensitivity = FALSE` disables marker-sensitivity weighting — recommended
    when scoring with learned weights to avoid double-penalizing shared markers.
- 🔧 Statistical soundness: removed input masking in `prepare_data_for_training()`
  (masked zeros let the model learn the mask pattern instead of expression).
- 🔧 Seurat v5 compatibility: use `Seurat::GetAssayData(x, layer = "data")`
  instead of direct `@assays$RNA@data` / `@assays$RNA@scale.data` access.
- 🧩 New suggested dependencies: `glmnet`, `ranger`, `xgboost` (required per
  method, optional overall).

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
