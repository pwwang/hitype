# Extract the per-type linear predictor coefficients of a glm fit

Aligned by name to c("(Intercept)", markers): an aliased term can be
absent from [`coef()`](https://rdrr.io/r/stats/coef.html), it then stays
0.

## Usage

``` r
model_coefs_lr(fit, markers, return_models)
```
