# Update sequencing rates in a surveillance design

Update sequencing rates in a surveillance design

## Usage

``` r
surv_update_rates(design, new_rates)
```

## Arguments

- design:

  A `surv_design` object.

- new_rates:

  Data frame with strata columns + `seq_rate`, or named numeric vector.

## Value

Updated `surv_design`.

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 4, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
new_r <- sim$population[c("region", "seq_rate")]
new_r$seq_rate <- new_r$seq_rate * 2
d2 <- surv_update_rates(d, new_r)
```
