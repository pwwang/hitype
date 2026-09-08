<!-- README.md is generated from README.Rmd. Please edit that file -->

# hitype <a href="https://pwwang.github.io/hitype/"><img src="man/figures/logo.png" align="right" height="139" alt="hitype website" /></a>

<!-- badges: start -->
<!-- badges: end -->

**Hi**erarchical and **hi**gh-resolution cell-type identification for
single-cell RNA-seq data inspired by
[ScType](https://github.com/IanevskiAleksandr/sc-type).

## Features

-   [x] Compatibility with
    [ScType](https://github.com/IanevskiAleksandr/sc-type)
-   [x] Hierarchical and high-resolution cell-type identification
-   [x] Learn marker weights from data with 7 methods: `uniform`,
    `correlation`, `lr`, `glmnet` (default), `rf`, `xgb`, `lrp`
-   [x] Cross-validated weight training (`cv_folds`) for stable,
    statistically sound weights
-   [x] Optional marker-sensitivity weighting and score normalization
    (`use_sensitivity`, `norm`)
-   [x] Speed optimization
-   [x] Support as an R package with unit tests

## Installation

You can install the development version of `hitype` like so:

``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("pwwang/hitype")
```

Optional packages for the individual weight-learning methods:
`glmnet` (method `glmnet`), `ranger` (method `rf`), `xgboost` (method
`xgb`), `keras` + `innsight` (method `lrp`).

## Quick start

### Prepare the dataset

See also
<https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#setup-the-seurat-object>

<details>
<summary>
Click to expand
</summary>

``` r
pbmc <- pbmc3k.SeuratData::pbmc3k
pbmc <- Seurat::UpdateSeuratObject(pbmc)
pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- Seurat::NormalizeData(pbmc)
pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- Seurat::ScaleData(pbmc, features = rownames(pbmc))
pbmc <- Seurat::RunPCA(pbmc, features = Seurat::VariableFeatures(object = pbmc))
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)
```

</details>

### Use as a Seurat extension

``` r
library(hitype)

markers <- data.frame(
    cellName = c(
        "Naive CD4+ T", "CD14+ Mono", "Memory CD4+", "B",
        "CD8+ T", "FCFR3A+ Mono", "NK", "DC", "Platelet"
    ),
    geneSymbolmore1 = c(
        "IL7R,CCR7",  "CD14,LYZ", "IL7R,S100A4", "MS4A1",
        "CD8A", "FCGR3A,MS4A7", "GNLY,NKG7", "FCER1A,CST3", "PPBP"
    ),
    geneSymbolmore2 = rep("", 9)
)

# Load gene sets
gs <- gs_prepare(markers)

# Assign cell types
obj <- RunHitype(pbmc, gs)

Seurat::DimPlot(obj, group.by = "hitype", label = TRUE, label.box = TRUE) +
  Seurat::NoLegend()
```

Compared to the manual marked cell types: <img
src="https://satijalab.org/seurat/articles/pbmc3k_tutorial_files/figure-html/labelplot-1.png"
style="width:75.0%" alt="Seurat manual marked cell types" />

See also
<https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#assigning-cell-type-identity-to-clusters>

### Use as standalone functions

``` r
scores <- hitype_score(Seurat::GetAssayData(pbmc, layer = "data"), gs)
cell_types <- hitype_assign(pbmc$seurat_clusters, scores, gs)
summary(cell_types)
```

You may see that we have exactly the same assignment in the Seurat
tutorial:

| Cluster ID | Markers       | Cell Type    |
|:-----------|:--------------|:-------------|
| 0          | IL7R, CCR7    | Naive CD4+ T |
| 1          | CD14, LYZ     | CD14+ Mono   |
| 2          | IL7R, S100A4  | Memory CD4+  |
| 3          | MS4A1         | B            |
| 4          | CD8A          | CD8+ T       |
| 5          | FCGR3A, MS4A7 | FCGR3A+ Mono |
| 6          | GNLY, NKG7    | NK           |
| 7          | FCER1A, CST3  | DC           |
| 8          | PPBP          | Platelet     |

`hitype_score()` accepts log-normalized data by default
(`scaled = FALSE`; it performs its own z-scoring) or pre-scaled data
with `scaled = TRUE`. When scoring with learned weights (see below),
use `norm = "weight", use_sensitivity = FALSE` — this normalizes
scores by the total marker weight and avoids double-penalizing shared
markers.

### Train marker weights

`train_weights()` learns per-marker weights from a labeled dataset
(clusters or cell types) with 7 backends:

| `method`     | Description                                        | Requires        |
|:-------------|:---------------------------------------------------|:----------------|
| `uniform`    | Equal weights (baseline)                           | —               |
| `correlation`| Correlation of each marker with the cluster label  | —               |
| `lr`         | Logistic regression (one-vs-rest) coefficients     | —               |
| `glmnet`     | Sparse elastic-net logistic regression (default)   | `glmnet`        |
| `rf`         | Random forest permutation importance               | `ranger`        |
| `xgb`        | XGBoost gain-based feature importance              | `xgboost`       |
| `lrp`        | Neural network + Layer-wise Relevance Propagation  | `keras`, `innsight` |

``` r
weights <- train_weights(
    path_to_gs = markers,
    exprs = pbmc,
    method = "glmnet",
    cv_folds = 5
)
gs <- gs_prepare(weights)
```

`cv_folds > 1` performs stratified cross-validation inside the training
data and averages the weights across folds for stability. **Do not
evaluate on data used for training** — train the weights on a train
split and score on a held-out split (or another dataset) to avoid
optimistic, circular results. See `vignette("train-marker-weights")`
for a full walkthrough including cross-dataset transfer.

### Find markers from your data

`find_markers()` discovers marker genes from a labeled dataset and
returns them directly in the database format consumed by `gs_prepare()`.
The default `method = "fc"` is dependency-light; `"seurat"` and
`"presto"` backends are also available:

``` r
markers <- find_markers(
    exprs = pbmc, clusters = pbmc$seurat_clusters, method = "fc"
)
weights <- train_weights(
    path_to_gs = markers, exprs = pbmc, method = "glmnet"
)
gs <- gs_prepare(weights)
```

**Fair-use warning:** markers and weights derived from the same dataset
are for exploratory use only. For publications, keep the marker database
fixed and follow the benchmark protocol (see `paper_plan.md`) to avoid
circular results.

## Documentation

<https://pwwang.github.io/hitype/>