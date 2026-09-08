# Find marker genes for cell types

Discovers marker genes for each cell type (cluster) in an expression
matrix or Seurat object, and returns them in the hitype marker database
format so the result can be passed directly to
[`gs_prepare()`](https://pwwang.github.io/hitype/reference/gs_prepare.md).
Three backends are available: the dependency-light fold-change method
(default), Seurat's `FindAllMarkers`, and presto's Wilcoxon rank-sum
test. For `method = "fc"`, the input expression matrix is expected to be
log-normalized.

## Usage

``` r
find_markers(
  exprs,
  clusters = NULL,
  method = c("fc", "seurat", "presto"),
  top = 20,
  min_log2fc = 0.25,
  min_pct = 0.1,
  only_pos = TRUE,
  include_negative = FALSE,
  level = 1
)
```

## Arguments

- exprs:

  A genes x cells expression matrix (plain matrix or
  dgCMatrix/dgTMatrix) or a Seurat object. For `method = "fc"`, the
  input is expected to be log-normalized.

- clusters:

  A named vector of cluster ids (names = cells) or a factor. If `exprs`
  is a Seurat object, defaults to `Seurat::Idents(exprs)`.

- method:

  The marker-finding method. One of:

  "fc"

  :   Fold-change based (default). No extra dependencies. Ranks genes by
      `log2fc * (pct_in - pct_out)`.

  "seurat"

  :   Seurat's `FindAllMarkers`. Requires a Seurat object.

  "presto"

  :   presto's Wilcoxon rank-sum test (`wilcoxauc`). Requires the presto
      package.

- top:

  Number of markers to return per cell type.

- min_log2fc:

  Minimum log2 fold change for a gene to be kept as a marker (when
  `only_pos = TRUE`).

- min_pct:

  Minimum fraction of cells in the cell type expressing the gene.

- only_pos:

  Only keep genes that are higher in the cell type than in the rest of
  the cells (positive markers). If `FALSE`, all genes passing `min_pct`
  are considered, ranked by score.

- include_negative:

  If `TRUE`, also fill the `geneSymbolmore2` column with the top
  down-regulated markers per cell type.

- level:

  The hierarchy level to write into the output database.

## Value

A data frame in the hitype db format with columns `cellName`,
`geneSymbolmore1`, `geneSymbolmore2` and `level`, directly consumable by
[`gs_prepare()`](https://pwwang.github.io/hitype/reference/gs_prepare.md).

## Examples

``` r
set.seed(1)
ngenes <- 40
ncells <- 60
exprs <- matrix(runif(ngenes * ncells, 0, 0.2), ngenes, ncells)
rownames(exprs) <- paste0("gene", seq_len(ngenes))
colnames(exprs) <- paste0("cell", seq_len(ncells))
clusters <- setNames(
    rep(c("Tcell", "Bcell", "Monocyte"), each = 20),
    colnames(exprs)
)
markers <- find_markers(exprs, clusters, method = "fc", top = 5)
head(markers)
#>   cellName                    geneSymbolmore1 geneSymbolmore2 level
#> 1    Tcell gene11,gene13,gene20,gene21,gene25                     1
#> 2    Bcell         gene7,gene14,gene17,gene38                     1
#> 3 Monocyte               gene21,gene26,gene38                     1
```
