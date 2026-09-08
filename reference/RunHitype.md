# Run hitype_assign for a Seurat object

Run hitype_assign for a Seurat object

## Usage

``` r
RunHitype(object, ...)

# Default S3 method
RunHitype(object, ...)

# S3 method for class 'Seurat'
RunHitype(
  object,
  gs = NULL,
  fallback = "Unknown",
  threshold = NULL,
  level_weights = function(l) 1/(10^(l - 1)),
  make_unique = FALSE,
  layer = "data",
  assay = NULL,
  scaled = FALSE,
  ...
)
```

## Arguments

- object:

  Seurat object

- ...:

  Additional arguments passed to the specific method.

- gs:

  The gene list prepared by
  [`gs_prepare`](https://pwwang.github.io/hitype/reference/gs_prepare.md)

- fallback:

  A fallback cell type if no cell type is assigned

- threshold:

  Confidence threshold as top1/top2 score ratio, passed to
  [`hitype_assign()`](https://pwwang.github.io/hitype/reference/hitype_assign.md).
  `NULL` (default) means no filtering.

- level_weights:

  The weights for each level of the hierarchy to calculate the final
  cell type score It should be either a numeric vector of length equal
  to the number of levels or a single numeric value to be used for all
  levels It can also be a function that takes the levels as input and
  returns a numeric vectors as the weights.

- make_unique:

  Whether to make the cell type names unique

- layer:

  The layer to use for `GetAssayData`

- assay:

  The assay to use for `GetAssayData`

- scaled:

  Whether the data from `GetAssayData` is scaled

## Value

The Seurat object with the cell types (named `hitype`) added to the
metadata
