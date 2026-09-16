# Prepare gene sets for hitype

Prepare gene sets for hitype

## Usage

``` r
gs_prepare(
  path_to_db_file,
  tissue_type = NULL,
  weight_encoding = function(x) x
)
```

## Arguments

- path_to_db_file:

  A data frame with markers or a path to a marker file. Two formats are
  supported and auto-detected:

  - The hitype/ScType db format (a data.frame or a tab-delimited text or
    excel file with the following columns):

    - `tissueType`: The tissue type of the cell type to be annotated.
      This column is required only if `tissue_type` is specified.

    - `cellName`: The name of the cell type to be annotated.

    - `nextLevels`: Possible next levels of the cell type to be
      annotated. Introduced by `hitype`, so that we can work with
      hierarchical cell. Levels are separated by `;`. Cell names at each
      level are separated by `,`. An exclamatory mark `!` at the
      beginning of a level means that the cell names at this level are
      mutually exclusive. If the levels are less than possible next
      levels, then the remaining levels are all possible next levels.
      See the example below.

    - `geneSymbolmore1`: The gene symbols of the marker genes that are
      expected to be expressed in the cell type to be annotated. The
      genes can be suffixed with one or more `+`. More `+` means higher
      expression level. For example, `CD3E++` means the gene `CD3E` is
      expected to be highly expressed in the cell type.

    - `geneSymbolmore2`: The gene symbols of the marker genes that are
      expected not to be expressed in the cell type to be annotated.

    - `level`: The levels of the cell names. Introduced by `hitype`, so
      that we can work with hierarchical cell names. Different levels of
      `cellName`s are predicted separately. For example, If we have
      `CD4` as level 1 and `Naive` as level 2, then our prediction for a
      cell type could be `CD4 Naive`. The levels should start from 1 and
      be consecutive.

  - The universal marker format (see
    <https://pwwang.github.io/biopipen/api/biopipen.ns.scrna/>), a long
    table with one row per gene per cell type:

    - `cell_type` (required): the cell type.

    - `gene` (required): the marker gene.

    - `direction`: `positive` or `negative` (aliases:
      `pos`/`neg`/`+`/`-`, case-insensitive). Optional. Defaults to
      `positive`.

    - `weight`: a numeric weight. Optional. Defaults to 1. A positive
      marker gets `abs(weight)`; a negative marker gets `-abs(weight)`.
      If no `direction` is given, the weight is used as-is (so signed
      weights are accepted).

    - `tissue`, `species`, `cancer`: optional, for filtering by
      `tissue_type`.

    - `level`: optional. The levels should start from 1 and be
      consecutive. Column names are case-insensitive and aliases are
      supported: `celltype`/`cellType`/`Type` for `cell_type`,
      `marker`/`gene_symbol` for `gene`, `sign` for `direction`, and
      `tissueType` for `tissue`. Text files with extensions
      `txt`/`tsv`/`csv`/`xlsx`/`xls` are read as tables;
      `rds`/`qs`/`qs2` files should contain a data.frame (`qs`/`qs2`
      require the qs2 package to be installed).

- tissue_type:

  The tissue type of the cell type to be annotated. For the db format,
  this requires the `tissueType` column in the marker gene database
  file; for the universal format, the `tissue` column (or its alias
  `tissueType`). If `tissue_type` is specified, then only the cell types
  in the specified tissue type will be used for annotation. If
  `tissue_type` is not specified, then all cell types in the marker gene
  database file will be used for annotation.

- weight_encoding:

  How to encoding the weights. By default, plus (+) for a positive
  weight; minus (-) for a negative weight and a star (\*) for zero
  weight. No sign indicates 1. One plus indicates 2, etc. You can use a
  function to encode these values. For example `function(x) x*x` to make
  the weights squared before they are involved in the computation.

## Value

A list with gene_sets and next_levels. The structure looks like:

      list(
         gene_sets = list(
             # level 1
             list( CD4 = list(markers = c(...), weights = c(...)), ... ),
             # level 2
             list( Naive = list(markers = c(...), weights = c(...)), ... )
         ),
         # All possible final cell names
         cell_names = list(CD4 = list(Naive = c("Activated", "Proliferating")))
      )

## Author

Matt Mulvahill, Panwen Wang

## Examples

``` r
# nextLevels example
# If we have the following cell types:
#   level  cellName  nextLevels
#   1      CD4       Naive,Memory
#   2      Naive
#   2      Memory
#   3      Activated
#   3      Proliferating
# Then possible final cell names are:
#   CD4 Naive Activated
#   CD4 Naive Proliferating
#   CD4 Memory Activated
#   CD4 Memory Proliferating
#
# If the `nextLevels` of CD4 is `!Naive`, then possible final cell names are:
#   CD4 Memory Activated
#   CD4 Memory Proliferating
#
# If the `nextLevels` of CD4 is `!`, then possible final cell names are:
#   CD4 Activated
#   CD4 Proliferating
#
# If the `nextLevels` of CD4 is `!Naive;!`, then possible final cell names
# are:
#   CD4 Memory

# A gene set in the universal marker format:
markers <- data.frame(
    cell_type = c("CD4", "CD4", "CD8", "CD8", "CD8"),
    gene = c("CD4", "IL7R", "CD8A", "CD8B", "GZMB"),
    direction = c("positive", "positive", "positive", "positive", "negative"),
    weight = c(2, 1, 3, 1, 1)
)
gs <- gs_prepare(markers)
gs$gene_sets[[1]]$CD4
#> $markers
#> [1] "CD4"  "IL7R"
#> 
#> $weights
#> [1] 2 1
#> 
gs$gene_sets[[1]]$CD8
#> $markers
#> [1] "CD8A" "CD8B" "GZMB"
#> 
#> $weights
#> [1]  3  1 -1
#> 
```
