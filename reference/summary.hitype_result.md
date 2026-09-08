# Summarize the hitype_result object

Summarize the hitype_result object

## Usage

``` r
# S3 method for class 'hitype_result'
summary(
  object,
  top = 1,
  level_weights = function(l) 1/(2^(l - 1)),
  make_unique = FALSE,
  ...
)
```

## Arguments

- object:

  A hitype_result object

- top:

  The top assigned cell types to return for each cluster

- level_weights:

  The weights for each level of the hierarchy to calculate the final
  cell type score It should be either a numeric vector of length equal
  to the number of levels or a single numeric value to be used for all
  levels It can also be a function that takes the levels as input and
  returns a numeric vectors as the weights.

- make_unique:

  Whether to make the cell type names unique

- ...:

  Additional arguments passed to the specific method.

## Value

The summary of the hitype_result object
