# survinger

<!-- badges: start -->
<!-- badges: end -->

**Design-adjusted inference for pathogen lineage surveillance.**

survinger provides tools for optimizing sequencing resource allocation and
estimating pathogen lineage prevalence under real-world genomic surveillance
conditions where sequencing rates are unequal across regions, sample sources
differ in cost and sensitivity, and reporting delays cause right-truncation
in recent data.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("CuiweiG/survinger")
```

## Quick example

```r
library(survinger)

# Simulate surveillance data (or use your own)
sim <- surv_simulate(n_regions = 5, n_weeks = 12, seed = 42)

# Create design object
design <- surv_design(
  data = sim$sequences,
  strata = ~ region,
  sequencing_rate = sim$population[c("region", "seq_rate")],
  population = sim$population
)

# Design-weighted prevalence (corrects for unequal sequencing)
weighted <- surv_lineage_prevalence(design, "BA.2.86")

# Naive prevalence (biased baseline)
naive <- surv_naive_prevalence(design, "BA.2.86")

# Compare
surv_compare_estimates(weighted, naive)

# Optimal allocation of 500 sequences
surv_optimize_allocation(design, "min_mse", total_capacity = 500)

# Delay correction
delay <- surv_estimate_delay(design)
adjusted <- surv_adjusted_prevalence(design, delay, "BA.2.86")
```

## Core functions

| Function | Purpose |
|----------|---------|
| `surv_design()` | Create surveillance design object |
| `surv_optimize_allocation()` | Optimize sequencing allocation |
| `surv_lineage_prevalence()` | Design-weighted prevalence (HT/Hajek/PS) |
| `surv_estimate_delay()` | Fit reporting delay distribution |
| `surv_nowcast_lineage()` | Nowcast recent incomplete counts |
| `surv_adjusted_prevalence()` | Combined design + delay correction |
| `surv_detection_probability()` | Variant detection power |
| `surv_compare_estimates()` | Weighted vs naive comparison plot |

## How it differs from phylosamp

**phylosamp** answers: *"How many sequences do I need?"*

**survinger** answers: *"How do I allocate my fixed capacity, and how do I correct the resulting estimates?"*

They are complementary, not competing.

## Citation

If you use survinger in published work, please cite:

> Gao, A. (2026). survinger: Design-Adjusted Inference for Pathogen
> Lineage Surveillance. R package version 0.1.0.
> https://github.com/CuiweiG/survinger

## License

MIT
