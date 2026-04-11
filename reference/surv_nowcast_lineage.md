# Nowcast lineage counts correcting for reporting delays

Nowcast lineage counts correcting for reporting delays

## Usage

``` r
# S3 method for class 'surv_nowcast'
print(x, ...)

# S3 method for class 'surv_nowcast'
as.data.frame(x, ...)

surv_nowcast_lineage(
  design,
  delay_fit,
  lineage = NULL,
  time = "epiweek",
  horizon = 4L,
  ref_date = NULL,
  method = c("direct", "em")
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

  Character or `NULL`. Target lineage.

- time:

  Character. Default `"epiweek"`.

- horizon:

  Integer. Recent periods to nowcast. Default 4.

- ref_date:

  Date or `NULL`. Default max report date.

- method:

  Character: `"direct"` (default) or `"em"`.

## Value

Invisibly returns the input object.

A `surv_nowcast` object.

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 12, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
fit <- surv_estimate_delay(d)
nc <- surv_nowcast_lineage(d, fit, "BA.2.86")
print(nc)
#> ── Delay-Adjusted Nowcast ──────────────────────────────────────────────────────
#> Method: "direct"
#> Truncation window: 4 periods
#> Lineage: "BA.2.86"
#> 
#> # A tibble: 9 × 10
#>   time     n_observed median_collection days_ago report_prob n_estimated      se
#>   <chr>         <int> <date>               <int>       <dbl>       <dbl>   <dbl>
#> 1 2024-W01          2 2024-01-06             107       1.000        2.00 2.05e-5
#> 2 2024-W02          1 2024-01-12             101       1.000        1.00 3.00e-5
#> 3 2024-W06          2 2024-02-07              74       1.000        2.00 1.07e-3
#> 4 2024-W07          2 2024-02-17              64       1.000        2.00 3.44e-3
#> 5 2024-W08          6 2024-02-22              59       1.000        6.00 1.06e-2
#> 6 2024-W09          6 2024-03-02              51       1.000        6.00 2.64e-2
#> 7 2024-W10         11 2024-03-06              47       1.000       11.0  5.59e-2
#> 8 2024-W11         14 2024-03-13              39       0.998       14.0  1.52e-1
#> 9 2024-W12         17 2024-03-21              32       0.993       17.1  3.54e-1
#> # ℹ 3 more variables: ci_lower <dbl>, ci_upper <dbl>, is_nowcast <lgl>
```
