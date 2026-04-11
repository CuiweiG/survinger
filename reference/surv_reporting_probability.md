# Compute cumulative reporting probability

Compute cumulative reporting probability

## Usage

``` r
surv_reporting_probability(delay_fit, delta, stratum = NULL)
```

## Arguments

- delay_fit:

  A `surv_delay_fit` object.

- delta:

  Integer vector. Days since collection.

- stratum:

  Character or `NULL`. Default `NULL`.

## Value

Numeric vector of probabilities.

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 12, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
fit <- surv_estimate_delay(d)
surv_reporting_probability(fit, delta = c(7, 14, 21))
#> [1] 0.4058020 0.7800219 0.9357006
```
