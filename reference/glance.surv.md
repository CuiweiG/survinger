# One-row summary of survinger model

One-row summary of survinger model

## Usage

``` r
# S3 method for class 'surv_prevalence'
glance(x, ...)

# S3 method for class 'surv_delay_fit'
glance(x, ...)

# S3 method for class 'surv_adjusted'
glance(x, ...)
```

## Arguments

- x:

  A survinger result object.

- ...:

  Additional arguments (currently unused).

## Value

A single-row tibble with key summary statistics.

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
prev <- surv_lineage_prevalence(d, "BA.2.86")
glance(prev)
#> # A tibble: 1 × 8
#>   lineage method n_periods mean_prevalence median_prevalence mean_se
#>   <chr>   <chr>      <int>           <dbl>             <dbl>   <dbl>
#> 1 BA.2.86 hajek         10          0.0403            0.0325  0.0220
#> # ℹ 2 more variables: mean_effective_n <dbl>, conf_level <dbl>
```
