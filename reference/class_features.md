# The features a cell type reports

The cell type's own markers when `class_markers` is given, every
candidate feature otherwise (the previous behaviour).

## Usage

``` r
class_features(class_markers, ct, markers)
```

## Arguments

- class_markers:

  A named list of the marker genes of each cell type, or `NULL`

- ct:

  The cell type

- markers:

  All candidate features

## Value

The feature vector of `ct`
