# Variant detection probability under current design

Variant detection probability under current design

## Usage

``` r
surv_detection_probability(
  design,
  true_prevalence,
  delay_fit = NULL,
  n_periods = 1L,
  detection_threshold = 1L
)
```

## Arguments

- design:

  A `surv_design` object.

- true_prevalence:

  Numeric in (0,1).

- delay_fit:

  Optional `surv_delay_fit`.

- n_periods:

  Integer. Accumulation periods. Default 1.

- detection_threshold:

  Integer. Min detections. Default 1.

## Value

A list with `overall`, `cumulative`, `by_stratum`, `parameters`.

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
surv_detection_probability(d, 0.01)
#> $overall
#> [1] 0.5336508
#> 
#> $cumulative
#> # A tibble: 1 × 2
#>   period p_detect
#>    <int>    <dbl>
#> 1      1    0.534
#> 
#> $by_stratum
#> # A tibble: 3 × 3
#>   stratum_id n_seq_per_period p_detect
#>        <int>            <dbl>    <dbl>
#> 1          1             22.4    0.202
#> 2          2             17.9    0.165
#> 3          3             35.6    0.301
#> 
#> $parameters
#> $parameters$true_prevalence
#> [1] 0.01
#> 
#> $parameters$n_periods
#> [1] 1
#> 
#> $parameters$detection_threshold
#> [1] 1
#> 
#> $parameters$delay_adjusted
#> [1] FALSE
#> 
#> 
```
