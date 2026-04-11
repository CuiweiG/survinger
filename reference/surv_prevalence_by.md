# Estimate prevalence by subgroup

Applies
[`surv_lineage_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md)
within subgroups defined by a grouping variable. Analogous to
[`survey::svyby()`](https://rdrr.io/pkg/survey/man/svyby.html) for
stratified survey analysis.

## Usage

``` r
surv_prevalence_by(
  design,
  lineage,
  by,
  time = "epiweek",
  method = "hajek",
  conf_level = 0.95
)
```

## Arguments

- design:

  A `surv_design` object.

- lineage:

  Character. Target lineage.

- by:

  Character. Column name to group by (e.g., `"region"` or
  `"source_type"`).

- time:

  Character. Time aggregation. Default `"epiweek"`.

- method:

  Character. Estimation method. Default `"hajek"`.

- conf_level:

  Numeric. Default 0.95.

## Value

A tibble with columns: group, time, lineage, prevalence, se, ci_lower,
ci_upper, n_obs, effective_n.

## See also

[`surv_lineage_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md),
[`surv_filter()`](https://cuiweig.github.io/survinger/reference/surv_filter.md)

## Examples

``` r
sim <- surv_simulate(n_regions = 4, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")],
                 sim$population, source_type = "source_type")
surv_prevalence_by(d, "BA.2.86", by = "region")
#> Warning: Lineage "BA.2.86" not found in data.
#> ℹ Available lineages: BA.5, Other, XBB.1.5
#> ℹ Did you mean "BA.5"?
#> # A tibble: 40 × 9
#>    group    time    lineage prevalence    se ci_lower ci_upper n_obs effective_n
#>    <chr>    <chr>   <chr>        <dbl> <dbl>    <dbl>    <dbl> <int>       <dbl>
#>  1 Region_A 2024-W… BA.2.86          0     0 2.78e-17    0.434     5           5
#>  2 Region_A 2024-W… BA.2.86          0     0 0           0.390     6           6
#>  3 Region_A 2024-W… BA.2.86          0     0 0           0.658     2           2
#>  4 Region_A 2024-W… BA.2.86          0     0 0           0.490     4           4
#>  5 Region_A 2024-W… BA.2.86          0     0 0           0.490     4           4
#>  6 Region_A 2024-W… BA.2.86          0     0 2.78e-17    0.434     5           5
#>  7 Region_A 2024-W… BA.2.86          0     0 0           0.658     2           2
#>  8 Region_A 2024-W… BA.2.86          0     0 2.78e-17    0.434     5           5
#>  9 Region_A 2024-W… BA.2.86          0     0 0           0.490     4           4
#> 10 Region_A 2024-W… BA.2.86          0     0 2.78e-17    0.434     5           5
#> # ℹ 30 more rows
```
