<h1 align="center">survinger</h1>
<p align="center"><em>Design-adjusted inference for pathogen lineage surveillance<br>under unequal sequencing and reporting delays</em></p>

<p align="center">
<a href="https://github.com/CuiweiG/survinger"><img src="https://img.shields.io/badge/R--CMD--check-passing-brightgreen" alt="R-CMD-check" /></a>
<a href="https://github.com/CuiweiG/survinger"><img src="https://img.shields.io/badge/version-0.1.0-blue" alt="version" /></a>
<a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="MIT" /></a>
</p>

---

## The problem

Genomic surveillance systems sequence unevenly. Denmark sequences 12% of cases; Romania sequences 0.3%. If you estimate lineage prevalence by counting sequences, the result is dominated by Denmark — regardless of what is actually circulating across Europe.

**On real ECDC data, this produces up to 14 percentage points of error:**

<p align="center">
<img src="man/figures/figA_necessity.png" width="90%" />
</p>

The red shaded area is the bias eliminated by design weighting. survinger corrects this using Horvitz-Thompson / Hajek estimators with Wilson score confidence intervals.

## The bias is structured, not random

<p align="center">
<img src="man/figures/fig3_bias_heatmap.png" width="90%" />
</p>

Each country's bias depends on its sequencing rate *and* its local prevalence, and both change over time. Poland (under-sequenced, high prevalence) is systematically underweighted by naive methods. A single correction factor cannot fix this — you need per-stratum, per-period weights.

## The correction works

<p align="center">
<img src="man/figures/fig6_benchmark.png" width="85%" />
</p>

In controlled simulation (50 replicates × 6 inequality levels), the Hajek estimator maintains 0.6–2.5 pp absolute bias while the naive estimator reaches 3.2–8.7 pp. The advantage holds across all levels of sequencing inequality.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("CuiweiG/survinger")
```

## Quick example

```r
library(survinger)

# Create design from surveillance data
design <- surv_design(
  data = sequences, strata = ~ region,
  sequencing_rate = population[c("region", "seq_rate")],
  population = population
)

# Corrected prevalence (one line)
surv_lineage_prevalence(design, "BA.2.86")

# Or even simpler:
surv_estimate(data, ~ region, rates, pop, "BA.2.86")

# Full pipeline
delay <- surv_estimate_delay(design)
surv_adjusted_prevalence(design, delay, "BA.2.86")

# How should I allocate 500 sequences?
surv_optimize_allocation(design, "min_mse", total_capacity = 500)

# Is my system powerful enough?
surv_detection_probability(design, true_prevalence = 0.01)

# One-page diagnostic
surv_report(design)
```

## Core functions

| Function | Purpose |
|----------|---------|
| `surv_design()` | Create design with inverse-probability weights |
| `surv_lineage_prevalence()` | Hajek / HT / post-stratified prevalence |
| `surv_optimize_allocation()` | Neyman allocation (3 objectives) |
| `surv_estimate_delay()` | Right-truncation-corrected delay fitting |
| `surv_adjusted_prevalence()` | Combined design + delay correction |
| `surv_detection_probability()` | Variant detection power |
| `surv_report()` | Surveillance system diagnostic |
| `surv_quality()` | One-row quality metrics |
| `tidy()` / `glance()` | Tidyverse integration |

## How it differs from existing tools

| | phylosamp | survey | epinowcast | **survinger** |
|---|---|---|---|---|
| Question | How many? | General surveys | Bayesian nowcast | **Allocate + correct + nowcast** |
| Genomic-specific | ✓ | ✗ | Partial | **✓** |
| Allocation | ✗ | ✗ | ✗ | **✓ (3 objectives)** |
| Delay correction | ✗ | ✗ | ✓ | **✓** |
| Requires Stan | ✗ | ✗ | ✓ | **✗** |
| CRAN-friendly | ✓ | ✓ | ✗ | **✓** |

## Validated on real data

- **ECDC**: 99,093 sequences, 5 EU countries, 40-fold inequality
- **COG-UK**: 65,166 individual sequences, 4 UK nations
- Cross-validated against `survey::svymean` (exact match)
- Wilson CI coverage: 93.4% (Brown et al. 2001 target: 93–95%)
- Delay MLE recovery: 0.5% error at n = 5,000

## Vignettes

- `vignette("survinger")` — Quick start
- `vignette("allocation-optimization")` — Resource allocation
- `vignette("delay-correction")` — Delay estimation and nowcasting
- `vignette("real-world-ecdc")` — ECDC case study

## Citation

```bibtex
@Manual{survinger2026,
  title = {survinger: Design-Adjusted Inference for Pathogen Lineage Surveillance},
  author = {Cuiwei Gao},
  year = {2026},
  note = {R package version 0.1.0},
  url = {https://github.com/CuiweiG/survinger}
}
```

## License

MIT © 2026 Cuiwei Gao
