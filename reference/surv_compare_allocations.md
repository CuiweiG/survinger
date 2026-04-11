# Compare multiple allocation strategies

Compare multiple allocation strategies

## Usage

``` r
surv_compare_allocations(
  design,
  strategies = c("equal", "proportional", "min_mse", "max_detection", "min_imbalance"),
  total_capacity,
  target_prevalence = 0.01,
  ...
)
```

## Arguments

- design:

  A `surv_design` object.

- strategies:

  Character vector. Default includes all built-in.

- total_capacity:

  Integer. Total sequences.

- target_prevalence:

  Numeric. For detection objective.

- ...:

  Passed to
  [`surv_optimize_allocation()`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md).

## Value

A tibble comparing strategies.

## Examples

``` r
sim <- surv_simulate(n_regions = 4, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
surv_compare_allocations(d, total_capacity = 200)
#> # A tibble: 5 × 4
#>   strategy      total_mse detection_prob imbalance
#>   <chr>             <dbl>          <dbl>     <dbl>
#> 1 equal           0.00149          0.866 0.0499   
#> 2 proportional    0.00124          0.866 0.0000103
#> 3 min_mse         0.00124          0.866 0.0000103
#> 4 max_detection   0.00158          0.866 0.0703   
#> 5 min_imbalance   0.00124          0.866 0.0000103
```
