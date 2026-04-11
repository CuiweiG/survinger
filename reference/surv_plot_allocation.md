# Plot allocation plan

Plot allocation plan

## Usage

``` r
surv_plot_allocation(allocation)
```

## Arguments

- allocation:

  A `surv_allocation` object.

## Value

A ggplot2 object.

## Examples

``` r
sim <- surv_simulate(n_regions = 5, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
a <- surv_optimize_allocation(d, "min_mse", total_capacity = 500)
surv_plot_allocation(a)

```
