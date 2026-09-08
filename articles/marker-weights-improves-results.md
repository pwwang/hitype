# Marker weights improve cell type assignment

``` r
library(hitype)
```

Unweighted markers:

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

Prepare `pbmc` data:

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

Use unweighted markers:

``` r
gs <- gs_prepare(markers)
scores <- hitype_score(Seurat::GetAssayData(pbmc, layer = "data"), gs)
assigned <- hitype_assign(pbmc$seurat_clusters, scores, gs, threshold = 0.01)
out <- summary(assigned)
out$ManualAssignment <- markers$cellName
out
#> # A tibble: 9 × 6
#>   Level Cluster CellType      Score  Margin ManualAssignment
#>   <int> <fct>   <chr>         <dbl>   <dbl> <chr>           
#> 1     1 0       Naive CD4+ T  0.721  0.685  Naive CD4+ T    
#> 2     1 1       CD14+ Mono    2.28   1.20   CD14+ Mono      
#> 3     1 2       Memory CD4+   0.599  0.0927 Memory CD4+     
#> 4     1 3       B             2.09   2.17   B               
#> 5     1 4       NK            1.43   0.179  CD8+ T          
#> 6     1 5       FCFR3A+ Mono  3.42   2.44   FCFR3A+ Mono    
#> 7     1 6       NK            3.91   2.88   NK              
#> 8     1 7       DC            5.55   4.59   DC              
#> 9     1 8       Platelet     11.5   11.3    Platelet
```

Use trained weights (glmnet, 5-fold CV). Note the recommended scoring
settings for learned weights:
`norm = "weight", use_sensitivity = FALSE`:

``` r
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

gs <- gs_prepare(weights)
scores <- hitype_score(Seurat::GetAssayData(pbmc, layer = "data"), gs,
                       norm = "weight", use_sensitivity = FALSE)
assigned <- hitype_assign(pbmc$seurat_clusters, scores, gs, threshold = 0.01)
out <- summary(assigned)
out$ManualAssignment <- markers$cellName
out
#> # A tibble: 9 × 6
#>   Level Cluster CellType      Score  Margin ManualAssignment
#>   <int> <fct>   <chr>         <dbl>   <dbl> <chr>           
#> 1     1 0       Naive CD4+ T  0.535  0.500  Naive CD4+ T    
#> 2     1 1       CD14+ Mono    1.61   0.847  CD14+ Mono      
#> 3     1 2       Memory CD4+   0.469  0.0656 Memory CD4+     
#> 4     1 3       B             2.09   2.17   B               
#> 5     1 4       CD8+ T        1.25   0.239  CD8+ T          
#> 6     1 5       FCFR3A+ Mono  2.42   1.72   FCFR3A+ Mono    
#> 7     1 6       NK            2.76   2.04   NK              
#> 8     1 7       DC            3.93   3.25   DC              
#> 9     1 8       Platelet     11.5   11.3    Platelet
```

Use unweighted markers on `ifnb` dataset:

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

gs <- gs_prepare(markers)

scores <- hitype_score(Seurat::GetAssayData(ifnb, layer = "data"), gs)

assigned <- hitype_assign(ifnb$seurat_clusters, scores, gs, threshold = 0.01)

summary(assigned)
```

Use trained weights on `ifnb` dataset (cross-dataset transfer — the
weights were learned on `pbmc`, never on `ifnb`):

``` r
gs <- gs_prepare(weights)

scores <- hitype_score(Seurat::GetAssayData(ifnb, layer = "data"), gs,
                       norm = "weight", use_sensitivity = FALSE)

assigned <- hitype_assign(ifnb$seurat_clusters, scores, gs, threshold = 0.01)

summary(assigned)
```
