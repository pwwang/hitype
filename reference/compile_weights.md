# Compile the weights

Compile the weights

## Usage

``` r
compile_weights(
  weights,
  gs,
  level,
  range = c(-5, -1, 1, 5),
  format = c("universal", "db"),
  drop_zero = TRUE,
  pos_only = FALSE
)
```

## Arguments

- weights:

  A data frame with the weights

- gs:

  The gene sets

- level:

  The level of the gene sets

- range:

  The quantization range for `format = "db"` only: the 4-element form
  `c(-low, -high, low, high)` (default `c(-5, -1, 1, 5)`) scales the
  positive weights into `[low, high]` and the negative ones into
  `[-high, -low]` before rounding them to integers (the db format can
  only encode integer weights via `+`/`-` repeats). The universal output
  is never rescaled — it keeps the raw trained weights.

- format:

  The format of the output data frame. One of `"universal"` (default) or
  `"db"`.

- drop_zero:

  Whether to drop markers whose trained weight is exactly zero from the
  output (default `TRUE`). Zero weights mean the marker carries no
  learned direction for the cell type; dropping them (before any db
  quantization) makes zero mean "no direction" — the marker is absent
  from the cell type's list. Pass `FALSE` to keep every candidate
  marker.

- pos_only:

  Whether to keep only the markers whose trained weight is positive in
  the output (default `FALSE`). A marker trained with a negative weight
  anti-correlates with the cell type and is reported with
  `direction = "negative"`; pass `TRUE` to drop those rows (as well as
  any exactly-zero ones, so `drop_zero` is then redundant) and keep only
  the markers overexpressed in the cell type.

## Value

A data frame with the compiled weights in the universal marker format
(default) or the db format (`format = "db"`), consumable by
[`gs_prepare()`](https://pwwang.github.io/hitype/reference/gs_prepare.md).
