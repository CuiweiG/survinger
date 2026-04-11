# Compute surveillance quality metrics

Returns a single-row tibble of design quality indicators suitable for
inclusion in manuscripts. Analogous to `broom::glance()` but for the
surveillance design itself.

## Usage

``` r
surv_quality(design, target_lineage = NULL, target_prevalence = 0.01)
```

## Arguments

- design:

  A `surv_design` object.

- target_lineage:

  Character or `NULL`. Default `NULL` auto-selects.

- target_prevalence:

  Numeric. For detection calculation. Default 0.01.

## Value

A single-row tibble with columns: n_obs, n_strata, gini, rate_ratio,
effective_n, deff, detection_prob, mean_bias.

## See also

[`surv_report()`](https://cuiweig.github.io/survinger/reference/surv_report.md),
[`surv_design()`](https://cuiweig.github.io/survinger/reference/surv_design.md)

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
surv_quality(d)
#> # A tibble: 1 × 12
#>   n_obs n_strata n_lineages  gini rate_min rate_max rate_ratio effective_n  deff
#>   <int>    <int>      <int> <dbl>    <dbl>    <dbl>      <dbl>       <dbl> <dbl>
#> 1   759        3          4 0.356   0.0395    0.220        5.6         573  1.33
#> # ℹ 3 more variables: detection_at_1pct <dbl>, mean_abs_bias_pp <dbl>,
#> #   target_lineage <chr>
```
