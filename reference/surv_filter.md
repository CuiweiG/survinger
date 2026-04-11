# Subset a surveillance design by filter criteria

Creates a new `surv_design` object containing only sequences matching
the specified filter criteria.

## Usage

``` r
surv_filter(design, ...)
```

## Arguments

- design:

  A `surv_design` object.

- ...:

  Filter expressions passed to
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html).

## Value

A new `surv_design` object with filtered data.

## See also

[`surv_design()`](https://cuiweig.github.io/survinger/reference/surv_design.md)

## Examples

``` r
sim <- surv_simulate(n_regions = 5, n_weeks = 12, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
d_sub <- surv_filter(d, region %in% c("Region_A", "Region_B"))
print(d_sub)
#> ── Genomic Surveillance Design ─────────────────────────────────────────────────
#> Observations: 626
#> Strata: 5 (region)
#> Date range: 2024-01-01 to 2024-03-24
#> Lineages: 4
#> Weight range: [2.003, 10.653]
```
