# Pipe-friendly surveillance analysis

Convenience wrapper that creates a design, estimates prevalence, and
optionally applies delay correction in a single pipe-friendly call.
Designed for rapid exploratory analysis in interactive sessions.

## Usage

``` r
surv_estimate(
  data,
  strata,
  sequencing_rate,
  population,
  lineage,
  correct_delay = FALSE,
  method = "hajek",
  ...
)
```

## Arguments

- data:

  Data frame of sequence records.

- strata:

  One-sided formula for stratification.

- sequencing_rate:

  Sequencing rate specification (see
  [`surv_design()`](https://cuiweig.github.io/survinger/reference/surv_design.md)).

- population:

  Population data frame.

- lineage:

  Character. Target lineage to estimate.

- correct_delay:

  Logical. Apply delay correction? Default `FALSE`.

- method:

  Character. Prevalence method. Default `"hajek"`.

- ...:

  Additional arguments passed to
  [`surv_design()`](https://cuiweig.github.io/survinger/reference/surv_design.md).

## Value

A `surv_prevalence` or `surv_adjusted` object.

## See also

[`surv_design()`](https://cuiweig.github.io/survinger/reference/surv_design.md),
[`surv_lineage_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md)

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
# One-liner analysis:
result <- surv_estimate(
  data = sim$sequences, strata = ~ region,
  sequencing_rate = sim$population[c("region", "seq_rate")],
  population = sim$population,
  lineage = "BA.2.86"
)
print(result)
#> ── Lineage Prevalence Estimate ─────────────────────────────────────────────────
#> Lineage: "BA.2.86"
#> Method: "hajek"
#> Confidence level: 0.95
#> Time periods: 10
#> 
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
