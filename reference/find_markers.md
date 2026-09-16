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
  top = c(10, 10),
  min_log2fc = 0.25,
  min_pct = 0.1,
  against = NULL,
  max_pct_out = 0.75,
  pos_only = TRUE,
  level = 1,
  format = c("universal", "db")
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

  Number of markers to return per cell type: a single number applied to
  both positive and negative markers, or a length-2 vector
  `c(n_positive, n_negative)` budgeting each direction separately
  (default `c(10, 10)`). The negative budget is only used when
  `pos_only = FALSE`.

- min_log2fc:

  Minimum log2 fold change for a gene to be kept as a marker (when
  `pos_only = TRUE`).

- min_pct:

  Minimum fraction of cells in the cell type expressing the gene.

- against:

  Restrict the reference group a cell type is compared against (by
  default one-vs-rest, i.e. all other cells). One of:

  `NULL`

  :   One-vs-rest (default).

  a character vector of cell types

  :   Only the cells of those types form the reference group: `pct_out`,
      `log2fc` and the score are computed against them instead of all
      the other cells, which recovers sibling-specific markers shared
      with the rest of the cells. A type listed in `against` cannot use
      its own cells as the reference, so they are excluded from its
      reference group; if no other listed type remains, no markers are
      returned for it.

  `"nearest"`

  :   For each cell type, the single most similar other type is used as
      its reference group. Similarity is the Pearson correlation between
      the mean expression profiles of the two types (computed once per
      call over the input matrix), and the resolved type is then used
      exactly like an explicit `against` entry.

  An error is raised if `against` names a cell type not present in the
  data (available types are listed) or if it names the only cell type in
  the data (a type cannot be compared against itself). Negative markers
  (`pos_only = FALSE`) are extracted as before as genes low in the type
  vs the rest of the cells, but when `against` is set they must ALSO
  satisfy `log2fc <= -min_log2fc` against the reference group, so that
  sibling-specific negative markers are found. For `method = "seurat"`,
  `against` is implemented with per-type two-group
  `Seurat::FindMarkers(ident.1 = type, ident.2 = against)` calls
  (`"nearest"` is resolved to the per-type reference types first);
  `method = "presto"` does not support `against` and raises an error.

- max_pct_out:

  Drop positive-marker candidates expressed in more than this fraction
  of ALL other cells (the one-vs-rest `pct_out`, computed over all cells
  not of the type even when `against` is set), guarding against
  pan-lineage genes up-regulated in most other cell types. Negative
  markers are not subject to this guard. A single aggregated warning
  reports the number of dropped candidates per cell type. Default
  `0.75`.

- pos_only:

  Only keep genes that are higher in the cell type than in the rest of
  the cells (positive markers). If `FALSE`, all genes passing `min_pct`
  are considered, ranked by score.

- level:

  The hierarchy level to write into the output data frame.

- format:

  The format of the output data frame. One of `"universal"` (default) or
  `"db"` (the hitype/ScType wide format).

## Value

A data frame in the universal marker format (default) with columns
`cell_type`, `gene`, `direction` and `level`, or in the hitype db format
(`format = "db"`) with columns `cellName`, `geneSymbolmore1`,
`geneSymbolmore2` and `level`. Both are directly consumable by
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
#> Warning: Dropped 16 positive-marker candidate(s) expressed in more than 75% of other cells (max_pct_out): Tcell: 9, Bcell: 4, Monocyte: 3
head(markers)
#> NULL
```
