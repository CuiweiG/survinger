# Compute design effect over time

Compute design effect over time

## Usage

``` r
surv_design_effect(weighted, naive)
```

## Arguments

- weighted:

  A `surv_prevalence` object.

- naive:

  A `surv_prevalence` object.

## Value

A tibble with time, deff, and bias_correction columns.

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
w <- surv_lineage_prevalence(d, "BA.2.86")
n <- surv_naive_prevalence(d, "BA.2.86")
surv_design_effect(w, n)
#> # A tibble: 10 × 7
#>    time     prev_w   se_w prev_n   se_n  deff bias_correction
#>    <chr>     <dbl>  <dbl>  <dbl>  <dbl> <dbl>           <dbl>
#>  1 2024-W01 0.0427 0.0318 0.0270 0.0189 2.85          0.0157 
#>  2 2024-W02 0.0366 0.0366 0.0172 0.0171 4.59          0.0194 
#>  3 2024-W03 0      0      0      0      0             0      
#>  4 2024-W04 0      0      0      0      0             0      
#>  5 2024-W05 0      0      0      0      0             0      
#>  6 2024-W06 0.0283 0.0244 0.0244 0.0170 2.05          0.00393
#>  7 2024-W07 0.0218 0.0152 0.0244 0.0170 0.797        -0.00259
#>  8 2024-W08 0.0605 0.0299 0.0769 0.0302 0.983        -0.0164 
#>  9 2024-W09 0.0907 0.0420 0.0870 0.0339 1.53          0.00376
#> 10 2024-W10 0.123  0.0403 0.128  0.0360 1.25         -0.00521
```
