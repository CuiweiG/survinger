# Compute naive (unweighted) lineage prevalence

Simple proportion without design correction. Useful as baseline.

## Usage

``` r
surv_naive_prevalence(design, lineage, time = "epiweek", conf_level = 0.95)
```

## Arguments

- design:

  A `surv_design` object.

- lineage:

  Character. Target lineage.

- time:

  Character. Time aggregation. Default `"epiweek"`.

- conf_level:

  Numeric. Default 0.95.

## Value

A `surv_prevalence` object with `method = "naive"`.

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
naive <- surv_naive_prevalence(d, "BA.2.86")
```
