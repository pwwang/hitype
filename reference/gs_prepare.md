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

  A data frame with markers or Path to the marker gene database file, it
  should be a tab-delimited text or excel file with the following
  columns:

  - `tissueType`: The tissue type of the cell type to be annotated. This
    column is required only if `cell_type` is specified.

  - `cellName`: The name of the cell type to be annotated.

  - `nextLevels`: Possible next levels of the cell type to be annotated.
    Introduced by `hitype`, so that we can work with hierarchical cell.
    Levels are separated by `;`. Cell names at each level are separated
    by `,`. An exclamatory mark `!` at the beginning of a level means
    that the cell names at this level are mutually exclusive. If the
    levels are less than possible next levels, then the remaining levels
    are all possible next levels. See the example below.

  - `geneSymbolmore1`: The gene symbols of the marker genes that are
    expected to be expressed in the cell type to be annotated. The genes
    can be suffixed with one or more `+`. More `+` means higher
    expression level. For example, `CD3E++` means the gene `CD3E` is
    expected to be highly expressed in the cell type.

  - `geneSymbolmore2`: The gene symbols of the marker genes that are
    expected not to be expressed in the cell type to be annotated.

  - `level`: The levels of the cell names. Introduced by `hitype`, so
    that we can work with hierarchical cell names. Different levels of
    `cellName`s are predicted separately. For example, If we have `CD4`
    as level 1 and `Naive` as level 2, then our prediction for a cell
    type could be `CD4 Naive`. The levels should start from 1 and be
    consecutive. \#'

- tissue_type:

  The tissue type of the cell type to be annotated. This requires the
  `tissueType` column in the marker gene database file. If `tissue_type`
  is specified, then only the cell types in the specified tissue type
  will be used for annotation. If `tissue_type` is not specified, then
  all cell types in the marker gene database file will be used for
  annotation.

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
```
