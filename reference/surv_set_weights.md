# Override design weights with custom values

Override design weights with custom values

## Usage

``` r
surv_set_weights(design, weights)
```

## Arguments

- design:

  A `surv_design` object.

- weights:

  Numeric vector of length equal to number of strata.

## Value

Updated `surv_design`.

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 4, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
d2 <- surv_set_weights(d, rep(1.0, d$n_strata))
```
