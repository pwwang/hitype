
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
#>   Cluster                      CellType
#> 0       0            CD4 Th17 Activated
#> 1       1           CD8 Naïve Activated
#> 2       2            CD4 Tscm Inhibited
#> 3       3       CD8 Tumor Recirculating
#> 4       4                CD4 Treg Naïve
#> 5       5  CD8 Tte Terminally Exhausted
#> 6       6                      MAIT Tem
#> 7       7 CD4 Naïve Precursor Exhausted
#> 8       8      MAIT Tumor Recirculating
```
