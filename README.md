
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hitype

<!-- badges: start -->
<!-- badges: end -->

**Hi**erarchical and **hi**gh-resolution cell-type identification for
single-cell RNA-seq data based on
[ScType](https://github.com/IanevskiAleksandr/sc-type).

## Features

-   [x] Compatibility with
    [ScType](https://github.com/IanevskiAleksandr/sc-type)
-   [x] Hierarchical and high-resolution cell-type identification
-   [x] Speed optimization
-   [x] Support as an R package with unit tests

## Installation

You can install the development version of hitype like so:

``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("pwwang/hitype")
```

## Quick start

``` r
library(hitype)

# Load gene sets
gs <- gs_prepare(hitypedb_tcell)

# Load expression data
pbmc3kt <- readRDS(url(
  "https://www.dropbox.com/scl/fi/pyizrlwuklt6g9yrgf51p/pbmc3kt.rds?rlkey=fz6t9qqjjf5n8dr08vv6rhyye&dl=1"
))

# Calculate cell type scores
scores <- suppressWarnings(  # Ignore non-exist genes
  hitype_score(pbmc3kt@assays$RNA@scale.data, gs, scaled = TRUE)
)

# Assign cell types
cell_types <- hitype_assign(pbmc3kt$seurat_clusters, scores, gs)
print(cell_types)
#>   Cluster CellType
#> 0       0  Unknown
#> 1       1  Unknown
#> 2       2  Unknown
#> 3       3  Unknown
#> 4       4  Unknown
#> 5       5  Unknown
#> 6       6  Unknown
#> 7       7  Unknown
#> 8       8  Unknown
```
