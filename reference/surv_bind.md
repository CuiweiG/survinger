# Combine multiple prevalence estimates

Bind results from different lineages or methods into a single tibble for
comparison plots and tables.

## Usage

``` r
surv_bind(...)
```

## Arguments

- ...:

  `surv_prevalence` objects to combine.

## Value

A tibble with a `source` column identifying each input.

## See also

[`surv_lineage_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md),
[`surv_sensitivity()`](https://cuiweig.github.io/survinger/reference/surv_sensitivity.md)

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
p1 <- surv_lineage_prevalence(d, "BA.5")
p2 <- surv_lineage_prevalence(d, "XBB.1.5")
combined <- surv_bind(p1, p2)
head(combined)
#> # A tibble: 6 × 11
#>   time    lineage n_obs n_target prevalence     se ci_lower ci_upper effective_n
#>   <chr>   <chr>   <int>    <int>      <dbl>  <dbl>    <dbl>    <dbl>       <dbl>
#> 1 2024-W… BA.5       74       41      0.574 0.0665    0.443    0.696        54.7
#> 2 2024-W… BA.5       58       31      0.546 0.0764    0.404    0.680        46.0
#> 3 2024-W… BA.5       74       48      0.626 0.0653    0.497    0.740        57.5
#> 4 2024-W… BA.5       78       38      0.478 0.0659    0.356    0.602        59.4
#> 5 2024-W… BA.5       78       45      0.571 0.0684    0.440    0.693        55.1
#> 6 2024-W… BA.5       82       41      0.437 0.0641    0.319    0.563        59.5
#> # ℹ 2 more variables: flag <chr>, source <chr>
```
