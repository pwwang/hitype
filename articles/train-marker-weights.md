# Train marker weights

You need to install the suggested packages for this vignette. The
default `method = "glmnet"` requires `glmnet`:

``` r
install.packages("glmnet")
devtools::install_github("pwwang/hitype")
```

Other backends have their own requirements (see
[`?train_weights`](https://pwwang.github.io/hitype/reference/train_weights.md)):
`ranger` (method `rf`), `xgboost` (method `xgb`), and `keras` +
`innsight` (method `lrp`). For `lrp` you would additionally run:

``` r
install.packages("keras")
install.packages("innsight")
keras::install_keras()
```

To train the weights for the markers, firstly, we have to have the
markers in the original `ScType` format:

``` r
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

markers
#>       cellName geneSymbolmore1 geneSymbolmore2
#> 1 Naive CD4+ T       IL7R,CCR7                
#> 2   CD14+ Mono        CD14,LYZ                
#> 3  Memory CD4+     IL7R,S100A4                
#> 4            B           MS4A1                
#> 5       CD8+ T            CD8A                
#> 6 FCFR3A+ Mono    FCGR3A,MS4A7                
#> 7           NK       GNLY,NKG7                
#> 8           DC     FCER1A,CST3                
#> 9     Platelet            PPBP
```

Prepare the data for training:

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
#> Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#> 
#> Number of nodes: 2638
#> Number of edges: 95927
#> 
#> Running Louvain algorithm...
#> Maximum modularity in 10 random starts: 0.8728
#> Number of communities: 9
#> Elapsed time: 0 seconds
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)
new_cluster_ids <- markers$cellName
names(new_cluster_ids) <- levels(pbmc)
pbmc <- Seurat::RenameIdents(pbmc, new_cluster_ids)
```

Train the weights with the recommended sparse logistic regression
(`glmnet`) and 5-fold cross-validation. `cv_folds > 1` stratifies the
training cells by cluster and averages the per-fold weights, making the
estimated weights stable and less prone to overfitting a single split:

``` r
library(hitype)

weights <- train_weights(
    path_to_gs = markers,
    exprs = pbmc,
    method = "glmnet",
    cv_folds = 5
)
#> # A tibble: 9 × 5
#>   Level Cluster      CellType     Score  Margin
#>   <dbl> <fct>        <chr>        <dbl>   <dbl>
#> 1     1 Naive CD4+ T Naive CD4+ T  2.71  2.58  
#> 2     1 CD14+ Mono   CD14+ Mono   10.2   5.39  
#> 3     1 Memory CD4+  Memory CD4+   2.79  0.922 
#> 4     1 B            B             7.80  8.14  
#> 5     1 CD8+ T       NK            5.60  0.0259
#> 6     1 FCFR3A+ Mono FCFR3A+ Mono 14.7  10.8   
#> 7     1 NK           NK           14.8   9.93  
#> 8     1 DC           DC           27.4  23.1   
#> 9     1 Platelet     Platelet     33.6  33.3
weights
#>       cellName level geneSymbolmore2    geneSymbolmore1
#> 1            B     1                           MS4A1+++
#> 2   CD14+ Mono     1                     CD14+++,LYZ+++
#> 3       CD8+ T     1                            CD8A+++
#> 4           DC     1                  FCER1A+++,CST3+++
#> 5 FCFR3A+ Mono     1                 FCGR3A+++,MS4A7+++
#> 6  Memory CD4+     1                  IL7R+++,S100A4+++
#> 7 Naive CD4+ T     1                    IL7R+++,CCR7+++
#> 8           NK     1                    GNLY+++,NKG7+++
#> 9     Platelet     1                            PPBP+++
```

Other methods are available via `method =` (e.g. `method = "lrp"` for
the deep-learning LRP backend, `method = "uniform"` for the equal-weight
baseline). For the `lr`, `glmnet`, and `lrp` backends, `cv_folds` is
used for cross-validation:

``` r
weights_lr <- train_weights(markers, pbmc, method = "lr", cv_folds = 5)
weights_lrp <- train_weights(markers, pbmc, method = "lrp", cv_folds = 5)
```

For a statistically sound evaluation, **do not score the trained weights
on the same cells used for training.** Either hold out a test split of
the training dataset, or transfer the weights to another dataset
(below).

Use the trained weights on `ifnb` dataset (cross-dataset transfer):

``` r
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))

ifnb <- ifnb.SeuratData::ifnb
ifnb <- Seurat::UpdateSeuratObject(ifnb)

# Keep only the cells that exist in the pbmc dataset
ifnb <- subset(ifnb, subset = seurat_annotations %in% c(
  "CD14 Mono",
  "DC",
  "CD4 Memory T",
  "CD4 Naive T",
  "CD8 T",
  "B",
  "NK"
))

# ifnb ships raw counts in SeuratData — normalize before scoring
ifnb <- Seurat::NormalizeData(ifnb)

ifnb@meta.data <- ifnb@meta.data %>%
  mutate(seurat_clusters = case_when(
    seurat_annotations == "CD14 Mono" ~ "CD14+ Mono",
    seurat_annotations == "CD4 Memory T" ~ "Memory CD4+",
    seurat_annotations == "CD4 Naive T" ~ "Naive CD4+ T",
    seurat_annotations == "CD8 T" ~ "CD8+ T",
    TRUE ~ seurat_annotations
  ))

gs <- gs_prepare(weights)

scores <- hitype_score(Seurat::GetAssayData(ifnb, layer = "data"), gs,
                       norm = "weight", use_sensitivity = FALSE)

assigned <- hitype_assign(ifnb$seurat_clusters, scores, gs, threshold = 0.01)

summary(assigned)
```

## Finding markers from your own data

[`find_markers()`](https://pwwang.github.io/hitype/reference/find_markers.md)
discovers marker genes from a labeled dataset (e.g.
`pbmc$seurat_clusters`) and returns them directly in the database format
consumed by
[`gs_prepare()`](https://pwwang.github.io/hitype/reference/gs_prepare.md),
with dependency-light fold-change (`fc`, default), `Seurat`, or `presto`
backends:

``` r
markers <- find_markers(pbmc, clusters = pbmc$seurat_clusters, method = "fc")
weights <- train_weights(path_to_gs = markers, exprs = pbmc, method = "glmnet")
gs <- gs_prepare(weights)
```

See
[`?find_markers`](https://pwwang.github.io/hitype/reference/find_markers.md)
for details.
