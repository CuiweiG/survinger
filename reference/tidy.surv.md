# Extract tidy estimates from survinger objects

Converts survinger result objects into tidy tibbles suitable for further
analysis with dplyr, ggplot2, or other tidyverse tools.

## Usage

``` r
# S3 method for class 'surv_prevalence'
tidy(x, ...)

# S3 method for class 'surv_nowcast'
tidy(x, ...)

# S3 method for class 'surv_adjusted'
tidy(x, ...)

# S3 method for class 'surv_allocation'
tidy(x, ...)

# S3 method for class 'surv_delay_fit'
tidy(x, ...)
```

## Arguments

- x:

  A survinger result object.

- ...:

  Additional arguments (currently unused).

## Value

A tibble.

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
prev <- surv_lineage_prevalence(d, "BA.2.86")
tidy(prev)
#> # A tibble: 10 × 10
#>    time   lineage n_obs n_target prevalence     se ci_lower ci_upper effective_n
#>    <chr>  <chr>   <int>    <int>      <dbl>  <dbl>    <dbl>    <dbl>       <dbl>
#>  1 2024-… BA.2.86    74        2     0.0427 0.0318  0.0128    0.133         54.7
#>  2 2024-… BA.2.86    58        1     0.0366 0.0366  0.00914   0.136         46.0
#>  3 2024-… BA.2.86    74        0     0      0       0         0.0627        57.5
#>  4 2024-… BA.2.86    78        0     0      0       0         0.0607        59.4
#>  5 2024-… BA.2.86    78        0     0      0       0         0.0651        55.1
#>  6 2024-… BA.2.86    82        2     0.0283 0.0244  0.00706   0.107         59.5
#>  7 2024-… BA.2.86    82        2     0.0218 0.0152  0.00478   0.0937        63.1
#>  8 2024-… BA.2.86    78        6     0.0605 0.0299  0.0226    0.152         59.2
#>  9 2024-… BA.2.86    69        6     0.0907 0.0420  0.0390    0.197         53.7
#> 10 2024-… BA.2.86    86       11     0.123  0.0403  0.0638    0.223         66.2
#> # ℹ 1 more variable: flag <chr>
```
