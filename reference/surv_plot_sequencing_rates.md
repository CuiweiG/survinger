# Plot sequencing rate inequality across strata

Plot sequencing rate inequality across strata

## Usage

``` r
surv_plot_sequencing_rates(design)
```

## Arguments

- design:

  A `surv_design` object.

## Value

A ggplot2 object.

## Examples

``` r
sim <- surv_simulate(n_regions = 5, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
surv_plot_sequencing_rates(d)

```
