# Estimate reporting delay distribution

Fits a parametric distribution to the delay between sample collection
and sequence reporting, accounting for right-truncation.

## Usage

``` r
# S3 method for class 'surv_delay_fit'
print(x, ...)

surv_estimate_delay(
  design,
  distribution = c("negbin", "poisson", "lognormal", "nonparametric"),
  strata = NULL,
  max_delay = 60L,
  ref_date = NULL
)
```

## Arguments

- x:

  Object to print.

- ...:

  Additional arguments (unused).

- design:

  A `surv_design` object with both collection and report dates.

- distribution:

  Character: `"negbin"` (default), `"poisson"`, `"lognormal"`, or
  `"nonparametric"`.

- strata:

  Optional one-sided formula for delay stratification.

- max_delay:

  Integer. Maximum plausible delay in days. Default 60.

- ref_date:

  Date. Reference date for right-truncation. Default is max report date
  in data.

## Value

Invisibly returns the input object.

A `surv_delay_fit` object.

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 12, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
fit <- surv_estimate_delay(d)
print(fit)
#> ── Reporting Delay Distribution ────────────────────────────────────────────────
#> Distribution: "negbin"
#> Strata: none (pooled)
#> Observations: 895
#> Mean delay: 10.1 days
#> 
#> # A tibble: 1 × 5
#>   stratum distribution    mu  size converged
#>   <chr>   <chr>        <dbl> <dbl> <lgl>    
#> 1 all     negbin        10.2  3.03 TRUE     
```
