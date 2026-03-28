<p align="center">
<strong style="font-size: 2em;">survinger</strong><br>
<em>Design-Adjusted Inference for Pathogen Lineage Surveillance</em>
</p>

<h1 align="center">survinger</h1>
<p align="center"><em>Design-adjusted inference for pathogen lineage surveillance<br>under unequal sequencing and reporting delays</em></p>

<!-- badges: start -->
<p align="center">
<a href="https://github.com/CuiweiG/survinger"><img src="https://img.shields.io/badge/R--CMD--check-passing-brightgreen" alt="R-CMD-check" /></a>
<a href="https://github.com/CuiweiG/survinger"><img src="https://img.shields.io/badge/version-0.1.0-blue" alt="version" /></a>
<a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="MIT" /></a>
<a href="https://github.com/CuiweiG/survinger"><img src="https://img.shields.io/badge/R-%E2%89%A5%204.1.0-blue" alt="R >= 4.1.0" /></a>
</p>
<!-- badges: end -->

---

## The problem

In pathogen genomic surveillance, **sequencing rates vary up to 40-fold across regions**. Naive lineage prevalence estimates — computed as simple sequence proportions — are dominated by high-sequencing regions regardless of true population-level trends. Reporting delays of 1–4 weeks further bias recent estimates downward through right-truncation.

**survinger** provides a unified statistical framework for bias correction.

---

## Why this package is necessary

### The cost of ignoring sequencing inequality
*Data: ECDC real sequences (5 EU countries, n = 99,093, 40-fold inequality)*

<p align="center">
<img src="man/figures/figA_necessity.png" width="90%" />
</p>

> **The red shaded area is the error you get from ignoring unequal sequencing.** On real ECDC data, the naive estimate (red) overestimates XBB.1.5-like prevalence by up to **14 percentage points** compared to the design-weighted estimate (blue). The bias is largest during the variant's peak — precisely when accurate estimates matter most. Mean absolute bias: 3.8 pp across all weeks.

### Where does the bias come from?
*Data: ECDC real sequences (5 EU countries, n = 99,093)*

<p align="center">
<img src="man/figures/figB_bias_source.png" width="80%" />
</p>

> **Mechanism:** Nations to the right of the 1x line are over-represented in sequencing relative to their population. If their local prevalence differs from the population mean, naive pooling produces biased estimates. survinger corrects this via inverse-probability weighting (Horvitz-Thompson / Hajek).

### Wilson intervals: valid coverage at any prevalence
*Data: Simulated emerging variant (~2% prevalence, small samples per week)*

<p align="center">
<img src="man/figures/figC_ci_comparison.png" width="85%" />
</p>

> **Method:** For rare or emerging variants, many weeks have zero observed sequences. Standard Wald CIs collapse to [0%, 0%] in these weeks (red dots, top panel — 6 of 15 weeks). survinger implements Wilson score intervals (bottom panel) which always provide nonzero width, maintaining valid 93.4% coverage. References: Wilson (1927) JASA; Agresti & Coull (1998).

---

## Validated on real surveillance data

Figures use **ECDC variant surveillance data** (5 EU countries, n = 99,093) and **COG-UK individual-level metadata** (4 UK nations, n = 65,166). Data sources: [ECDC Open Data](https://opendata.ecdc.europa.eu/covid19/virusvariant/), [COG-UK CLIMB](https://cog-uk.s3.climb.ac.uk/).

### Figure 1 · Sequencing inequality across countries
*Data: ECDC real data (5 EU countries, n = 99,093)*

<p align="center">
<img src="man/figures/fig1_inequality.png" width="85%" />
</p>

> **Observation:** Romania sequences 0.3% of cases vs Denmark at 12% — a **40-fold ratio**. This level of inequality is the norm in global genomic surveillance, and it makes naive prevalence estimates systematically biased.

### Figure 2 · Design weighting corrects systematic bias
*Data: ECDC real sequences (5 EU countries, n = 99,093)*

<p align="center">
<img src="man/figures/fig2_compare.png" width="85%" />
</p>

> **Result:** With 40-fold sequencing inequality (Gini = 0.54), the naive estimate diverges from the design-weighted estimate by **up to 14 percentage points** (mean 3.8 pp). The red shaded area is the bias eliminated by design weighting. The bias is largest during the variant's peak — precisely when accurate estimates matter most.

### Figure 3 · Bias structure varies by country and time
*Data: ECDC real sequences (5 EU countries, n = 99,093)*

<p align="center">
<img src="man/figures/fig3_bias_heatmap.png" width="90%" />
</p>

> **Interpretation:** Bias ranges from -37 pp to +45 pp across countries and weeks. Poland and Romania (lowest sequencing) show strong red (naive overestimates) while Denmark (highest sequencing) shows blue (naive underestimates). The pattern is time-varying and country-specific — no single correction factor suffices.

### Figure 4 · Delay-adjusted nowcasting (simulated data)

<p align="center">
<img src="man/figures/fig4_nowcast.png" width="85%" />
</p>

> **Method:** Right-truncation-corrected NegBin delay model. Demonstrated on simulated data (COG-UK does not publish upload dates). Recent weeks (▲) inflated by 1/F(Δ) where F is the estimated delay CDF.

### Figure 5 · Resource allocation optimization

<p align="center">
<img src="man/figures/fig5_allocation.png" width="85%" />
</p>

> **Finding:** MSE-optimal allocation concentrates resources in England (largest population, highest variance contribution). Equal allocation wastes 75% of NI's capacity given its small population.

### Figure 6 · Simulation benchmark: bias vs sequencing inequality
*Data: Controlled simulation (surv_simulate), 50 replicates × 6 Gini levels*

<p align="center">
<img src="man/figures/fig6_benchmark.png" width="85%" />
</p>

> **Key result:** Under heterogeneous prevalence (5%–30% across strata), the Hajek estimator maintains 0.6–2.5 pp absolute bias while the naive estimator reaches 3.2–8.7 pp. Both increase monotonically with inequality, but Hajek remains 3–8× lower. 50 replicates per Gini level; shaded bands show 95% CI.

### Figure 7 · Detection probability curve

<p align="center">
<img src="man/figures/fig7_detection.png" width="85%" />
</p>

> **Practical use:** `surv_detection_probability()` computes the probability of detecting ≥1 sequence of a variant at a given prevalence. With COG-UK sequencing volumes (n = 65,166 over 26 weeks), 50% weekly detection at 0.03%, 80% at 0.07%, 95% at 0.15% prevalence.

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
- `vignette("real-world-ecdc")` — Real-world case study (COG-UK / ECDC data)

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
