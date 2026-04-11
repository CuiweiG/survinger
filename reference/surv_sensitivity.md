# Sensitivity analysis across methods

Runs all three prevalence estimators and delay/no-delay variants on the
same design, producing a comparison table. Essential for robustness
checks in publications.

## Usage

``` r
surv_sensitivity(
  design,
  lineage,
  delay_fit = NULL,
  time = "epiweek",
  conf_level = 0.95
)
```

## Arguments

- design:

  A `surv_design` object.

- lineage:

  Character. Target lineage.

- delay_fit:

  Optional `surv_delay_fit` object. If provided, includes delay-adjusted
  estimates.

- time:

  Character. Default `"epiweek"`.

- conf_level:

  Numeric. Default 0.95.

## Value

A tibble with one row per method-time combination, columns: method,
time, prevalence, se, ci_lower, ci_upper.

## See also

[`surv_lineage_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md),
[`surv_adjusted_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_adjusted_prevalence.md)

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
surv_sensitivity(d, "BA.2.86")
#> # A tibble: 40 × 8
#>    method time     prevalence     se ci_lower ci_upper n_obs effective_n
#>    <fct>  <chr>         <dbl>  <dbl>    <dbl>    <dbl> <int>       <dbl>
#>  1 naive  2024-W01     0.0270 0.0189  0.00744   0.0933    74          74
#>  2 naive  2024-W02     0.0172 0.0171  0.00305   0.0914    58          58
#>  3 naive  2024-W03     0      0       0         0.0493    74          74
#>  4 naive  2024-W04     0      0       0         0.0469    78          78
#>  5 naive  2024-W05     0      0       0         0.0469    78          78
#>  6 naive  2024-W06     0.0244 0.0170  0.00671   0.0846    82          82
#>  7 naive  2024-W07     0.0244 0.0170  0.00671   0.0846    82          82
#>  8 naive  2024-W08     0.0769 0.0302  0.0357    0.158     78          78
#>  9 naive  2024-W09     0.0870 0.0339  0.0405    0.177     69          69
#> 10 naive  2024-W10     0.128  0.0360  0.0729    0.215     86          86
#> # ℹ 30 more rows
```
