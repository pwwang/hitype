# A frame is in the universal marker format iff it has (case-insensitively) a cell_type-family and a gene-family column. db frames never do.

A frame is in the universal marker format iff it has
(case-insensitively) a cell_type-family and a gene-family column. db
frames never do.

## Usage

``` r
is_universal_df(df)
```

## Arguments

- df:

  A data.frame of markers.

## Value

`TRUE` if the frame is in the universal marker format.
