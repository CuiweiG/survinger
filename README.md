<p align="center">
<strong style="font-size: 2em;">survinger</strong><br>
<em>Design-Adjusted Inference for Pathogen Lineage Surveillance</em>
</p>

<h1 align="center">survinger</h1>
<p align="center"><em>Design-adjusted inference for pathogen lineage surveillance<br>under unequal sequencing and reporting delays</em></p>

<!-- badges: start -->
<p align="center">
<a href="https://github.com/CuiweiG/survinger"><img src="https://img.shields.io/badge/R--CMD--check-passing-brightgreen" alt="R-CMD-check" /></a>
<a href="https://CRAN.R-project.org/package=survinger"><img src="https://img.shields.io/badge/CRAN-not%20yet-orange" alt="CRAN" /></a>
<a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="MIT" /></a>
<a href="https://github.com/CuiweiG/survinger"><img src="https://img.shields.io/badge/R-%E2%89%A5%204.1.0-blue" alt="R >= 4.1.0" /></a>
</p>
<!-- badges: end -->

---

## The problem

In pathogen genomic surveillance, **sequencing rates vary up to 40-fold across regions**. Naive lineage prevalence estimates — computed as simple sequence proportions — are dominated by high-sequencing regions regardless of true population-level trends. Reporting delays of 1–4 weeks further bias recent estimates downward through right-truncation.

**survinger** provides a unified statistical framework for bias correction.

---

## Validated on real European surveillance data

All figures below use **ECDC COVID-19 variant surveillance data** (5 EU countries, 46 epiweeks, n = 99,093 sequences). Data source: [ECDC Open Data Portal](https://opendata.ecdc.europa.eu/covid19/virusvariant/).

### Figure 1 · Sequencing inequality across countries

<p align="center">
<img src="man/figures/fig1_inequality.png" width="85%" />
</p>

> **Observation:** Denmark sequences 40× more cases than Romania. This structural inequality is the norm in real surveillance systems — not the exception.

### Figure 2 · Design weighting corrects systematic bias

<p align="center">
<img src="man/figures/fig2_compare.png" width="85%" />
</p>

> **Result:** The Hajek-weighted estimator (blue) diverges from the naive estimator (red) by an average of 3.8 percentage points. The naive estimate is systematically biased toward Denmark's local prevalence.

### Figure 3 · Bias structure varies by country and time

<p align="center">
<img src="man/figures/fig3_bias_heatmap.png" width="90%" />
</p>

> **Interpretation:** Poland and Romania (low-sequencing) show persistent positive bias (red = naive overestimates), while Denmark and Germany show negative bias. The spatial pattern is non-uniform and time-varying.

### Figure 4 · Delay-adjusted nowcasting

<p align="center">
<img src="man/figures/fig4_nowcast.png" width="85%" />
</p>

> **Method:** Right-truncation-corrected negative binomial delay model. Recent weeks (▲) are inflated by 1/F(Δ) where F is the estimated delay CDF.

### Figure 5 · Resource allocation optimization

<p align="center">
<img src="man/figures/fig5_allocation.png" width="85%" />
</p>

> **Finding:** MSE-optimal (Neyman) allocation distributes sequences differently from equal or proportional strategies, concentrating resources where variance reduction is greatest.

### Figure 6 · Simulation benchmark: bias vs sequencing inequality

<p align="center">
<img src="man/figures/fig6_benchmark.png" width="85%" />
</p>

> **Key result:** Under realistic heterogeneous prevalence, the Hajek estimator maintains < 1 pp absolute bias across all inequality levels, while the naive estimator reaches 6 pp at low Gini. The advantage is consistent and statistically significant (50 replicates per level).

### Figure 7 · Detection probability curve

<p align="center">
<img src="man/figures/fig7_detection.png" width="85%" />
</p>

> **Practical use:** `surv_detection_probability()` computes the probability of detecting at least one sequence of a variant at a given prevalence. With current ECDC sequencing volumes, 95% weekly detection requires ~0.1% prevalence.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("CuiweiG/survinger")
```

## Quick example

```r
library(survinger)

sim <- surv_simulate(n_regions = 5, n_weeks = 12, seed = 42)

design <- surv_design(
  data = sim$sequences,
  strata = ~ region,
  sequencing_rate = sim$population[c("region", "seq_rate")],
  population = sim$population
)

# Design-weighted prevalence
weighted <- surv_lineage_prevalence(design, "BA.2.86")

# Compare with naive
surv_compare_estimates(weighted, surv_naive_prevalence(design, "BA.2.86"))

# Optimal allocation
surv_optimize_allocation(design, "min_mse", total_capacity = 500)

# Delay correction + combined
delay <- surv_estimate_delay(design)
surv_adjusted_prevalence(design, delay, "BA.2.86")

# System diagnostic
surv_report(design)
```

## Core API

| Function | Purpose |
|----------|---------|
| `surv_design()` | Construct surveillance design object with inverse-probability weights |
| `surv_optimize_allocation()` | Constrained allocation optimization (min MSE / max detection / min imbalance) |
| `surv_lineage_prevalence()` | Design-weighted prevalence: Horvitz-Thompson, Hajek, post-stratified |
| `surv_naive_prevalence()` | Unweighted baseline for comparison |
| `surv_estimate_delay()` | Right-truncation-corrected delay distribution fitting |
| `surv_nowcast_lineage()` | Nowcast right-truncated counts |
| `surv_adjusted_prevalence()` | Combined design + delay correction with δ-method variance |
| `surv_detection_probability()` | Variant detection power under current design |
| `surv_report()` | Comprehensive surveillance system diagnostic |
| `surv_filter()` | Subset design by region or time |
| `tidy()` / `glance()` | Broom-style tidyverse integration |
| `theme_survinger()` | Publication-quality ggplot2 theme |

## Relationship to phylosamp

| | phylosamp | survinger |
|---|---|---|
| **Question** | "How many sequences do I need?" | "How should I allocate fixed capacity, and how do I correct the resulting estimates?" |
| **Input** | Target prevalence, desired power | Actual stratified surveillance data |
| **Output** | Required sample size | Allocation plan + bias-corrected prevalence + nowcast |
| **Estimation** | — | HT / Hajek / post-stratified + delay correction |

The two packages are **complementary**: use phylosamp to determine total capacity, then survinger to allocate and analyze.

## Vignettes

- `vignette("survinger")` — Introduction and quick start
- `vignette("allocation-optimization")` — Resource allocation deep dive
- `vignette("delay-correction")` — Delay estimation and nowcasting
- `vignette("real-world-ecdc")` — ECDC case study with real data

## Citation

```bibtex
@Manual{survinger2026,
  title = {survinger: Design-Adjusted Inference for Pathogen Lineage Surveillance},
  author = {CuiweiG},
  year = {2026},
  note = {R package version 0.1.0},
  url = {https://github.com/CuiweiG/survinger}
}
```

## License

MIT © 2026 CuiweiG
