# Generate a comprehensive surveillance system report

Produces a summary of the current surveillance design's strengths,
weaknesses, and recommendations.

## Usage

``` r
surv_report(design, target_lineage = NULL, target_prevalence = 0.01)
```

## Arguments

- design:

  A `surv_design` object.

- target_lineage:

  Character or `NULL`. Lineage to focus on. If `NULL`, uses the most
  common non-"Other" lineage.

- target_prevalence:

  Numeric. Assumed prevalence for detection calculations. Default 0.01.

## Value

Invisibly returns a named list of computed metrics including n_obs,
n_strata, rate_range, gini, effective_n, detection_prob, and mean_bias.

## See also

[`surv_design()`](https://cuiweig.github.io/survinger/reference/surv_design.md),
[`surv_lineage_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md),
[`surv_detection_probability()`](https://cuiweig.github.io/survinger/reference/surv_detection_probability.md),
[`surv_optimize_allocation()`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md)

## Examples

``` r
sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
surv_report(d)
#> 
#> ── Surveillance System Report ──────────────────────────────────────────────────
#> 
#> ── Design Structure ──
#> 
#> Total sequences: 759
#> Strata: 3
#> Period: 2024-01-01 to 2024-03-10
#> 
#> ── Sequencing Inequality ──
#> 
#> Rate range: 3.95% to 22.03%
#> Rate ratio: 5.6x
#> Gini coefficient: 0.356
#> ℹ Moderate inequality. Weighting recommended.
#> 
#> ── Estimation Impact ──
#> 
#> Effective sample size: 573 (75.4% of total)
#> Mean |weighted - naive| bias: 2.28 percentage points
#> 
#> ── Detection Power ──
#> 
#> Target: detect BA.5 at 1% prevalence
#> Weekly detection probability: 53.4%
#> Sequences needed for 95% detection: 299
```
