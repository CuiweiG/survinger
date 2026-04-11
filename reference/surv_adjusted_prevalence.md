# Combined design-weighted and delay-adjusted prevalence

Simultaneously corrects for unequal sequencing rates and
right-truncation from reporting delays.

## Usage

``` r
# S3 method for class 'surv_adjusted'
print(x, ...)

# S3 method for class 'surv_adjusted'
as.data.frame(x, ...)

surv_adjusted_prevalence(
  design,
  delay_fit,
  lineage,
  time = "epiweek",
  prevalence_method = "hajek",
  nowcast_method = "direct",
  conf_level = 0.95,
  bootstrap_n = 0L
)
```

## Arguments

- x:

  Object to print.

- ...:

  Additional arguments (unused).

- design:

  A `surv_design` object.

- delay_fit:

  A `surv_delay_fit` object.

- lineage:

  Character. Target lineage.

- time:

  Character. Default `"epiweek"`.

- prevalence_method:

  Character. Default `"hajek"`.

- nowcast_method:

  Character. Default `"direct"`.

- conf_level:

  Numeric. Default 0.95.

- bootstrap_n:

  Integer. 0 for delta method, \>0 for bootstrap. Default 0.

## Value

Invisibly returns the input object.

A `surv_adjusted` object.

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 12, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
delay <- surv_estimate_delay(d)
adj <- surv_adjusted_prevalence(d, delay, "BA.2.86")
print(adj)
#> ── Design-Weighted Delay-Adjusted Prevalence ───────────────────────────────────
#> Correction: "design:hajek+delay:direct"
#> 
#> # A tibble: 12 × 9
#>    time     lineage n_obs_raw n_obs_adjusted prevalence     se ci_lower ci_upper
#>    <chr>    <chr>       <int>          <dbl>      <dbl>  <dbl>    <dbl>    <dbl>
#>  1 2024-W01 BA.2.86        74           74.0     0.0432 0.0277   0        0.0975
#>  2 2024-W02 BA.2.86        58           58.0     0.0378 0.0283   0        0.0934
#>  3 2024-W03 BA.2.86        74           74.0     0      0        0        0     
#>  4 2024-W04 BA.2.86        78           78.0     0      0        0        0     
#>  5 2024-W05 BA.2.86        78           78.0     0      0        0        0     
#>  6 2024-W06 BA.2.86        82           82.0     0.0288 0.0219   0        0.0717
#>  7 2024-W07 BA.2.86        82           82.0     0.0213 0.0183   0        0.0572
#>  8 2024-W08 BA.2.86        78           78.0     0.0604 0.0312   0        0.122 
#>  9 2024-W09 BA.2.86        69           69.0     0.0914 0.0397   0.0137   0.169 
#> 10 2024-W10 BA.2.86        86           86.0     0.122  0.0406   0.0427   0.202 
#> # ℹ 2 more rows
#> # ℹ 1 more variable: mean_report_prob <dbl>
```
