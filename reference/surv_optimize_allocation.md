# Optimize sequencing allocation across strata

Given fixed total sequencing capacity, finds the optimal allocation
across strata that minimizes a specified objective function.

## Usage

``` r
# S3 method for class 'surv_allocation'
print(x, ...)

# S3 method for class 'surv_allocation'
as.data.frame(x, ...)

surv_optimize_allocation(
  design,
  objective = c("min_mse", "max_detection", "min_imbalance"),
  total_capacity,
  budget = NULL,
  min_per_stratum = 2L,
  target_lineage = NULL,
  target_prevalence = 0.01,
  cost_col = NULL
)
```

## Arguments

- x:

  Object to print.

- ...:

  Additional arguments (unused).

- design:

  A `surv_design` object.

- objective:

  Character. One of `"min_mse"`, `"max_detection"`, or
  `"min_imbalance"`.

- total_capacity:

  Integer. Total sequences available.

- budget:

  Numeric or `NULL`. Optional budget constraint.

- min_per_stratum:

  Integer. Minimum per stratum. Default 2.

- target_lineage:

  Character. Required for `"max_detection"`.

- target_prevalence:

  Numeric. Assumed prevalence for detection. Default 0.01.

- cost_col:

  Character or `NULL`. Column name for per-sequence cost.

## Value

Invisibly returns the input object.

A `surv_allocation` object.

## Examples

``` r
sim <- surv_simulate(n_regions = 4, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
a <- surv_optimize_allocation(d, "min_mse", total_capacity = 500)
print(a)
#> ── Optimal Sequencing Allocation ───────────────────────────────────────────────
#> Objective: min_mse
#> Total capacity: 500 sequences
#> Strata: 4
#> 
#> # A tibble: 4 × 3
#>   region   n_allocated proportion
#>   <chr>          <int>      <dbl>
#> 1 Region_A          44      0.088
#> 2 Region_B         178      0.356
#> 3 Region_C         176      0.352
#> 4 Region_D         102      0.204
```
