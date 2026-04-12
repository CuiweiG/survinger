# Summary

`survinger` is an R package that applies survey sampling methodology to
the problem of pathogen lineage prevalence estimation under unequal
genomic sequencing. It implements Horvitz–Thompson, Hajek, and
post-stratified estimators with Wilson score confidence intervals;
maximum likelihood estimation of right-truncated reporting delay
distributions for nowcasting; and Neyman-optimal resource allocation for
sequencing budgets. The package is validated against the `survey` R
package and against real European Centre for Disease Prevention and
Control (ECDC) surveillance data comprising 99,093 SARS-CoV-2 sequences
from five European countries.

# Statement of Need

Genomic surveillance for pathogen variants is now a cornerstone of
pandemic preparedness, yet the data it produces are systematically
biased by design. Sequencing rates vary enormously across jurisdictions
— Denmark sequences roughly 12% of positive COVID-19 cases while Romania
sequences 0.3%, a 40-fold difference \[@ecdc2023surveillance\]. When
lineage prevalence is estimated by naively counting sequences,
high-sequencing jurisdictions dominate the estimate regardless of their
population size. On real ECDC data, this produces estimates that deviate
from design-corrected values by up to 14 percentage points — a
difference large enough to alter variant risk assessments and trigger
inappropriate public health responses.

The statistical remedy is well-established. Survey sampling theory,
developed over decades for population surveys and census estimation,
provides a rigorous framework for inference from data collected with
unequal inclusion probabilities \[@horvitz1952generalization;
@cochran1977sampling\]. The Horvitz–Thompson estimator corrects for
unequal weights, the Hajek estimator improves efficiency through ratio
normalisation, and Neyman allocation minimises variance under fixed
total sample size constraints \[@neyman1934two\].

These methods are implemented in the `survey` R package
\[@lumley2004analysis\], which is designed for general-purpose survey
analysis. However, translating its interface to the specific data
structures, terminology, and operational concerns of genomic
surveillance — lineages rather than survey items, sequencing rates
rather than sampling fractions, reporting delays rather than
non-response — requires domain expertise that creates a barrier for
surveillance practitioners.

`survinger` bridges this gap. It implements the same statistical methods
in an interface expressed entirely in surveillance-native terminology,
adds right-truncation-corrected delay estimation and nowcasting
capabilities that are absent from general survey tools, and provides
sequencing resource allocation optimisation directly relevant to
surveillance programme management.

# Key Features

## Design-Weighted Prevalence Estimation

[`surv_design()`](https://cuiweig.github.io/survinger/reference/surv_design.md)
encodes the surveillance sampling design — strata definitions,
per-stratum sequencing rates, and population denominators — and computes
inverse-probability weights.
[`surv_lineage_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md)
estimates time-stratified prevalence for a target lineage using one of
three estimators:

- **Horvitz–Thompson**: Classical inverse-probability weighting
  \[@horvitz1952generalization\].
- **Hajek**: Ratio estimator normalising weights to sum to the known
  population total, producing more efficient estimates than
  Horvitz–Thompson when weights are variable.
- **Post-stratified**: Incorporates known population stratum sizes to
  improve precision.

All estimators produce Wilson score confidence intervals
\[@wilson1927probable\], which maintain proper coverage at extreme
proportions (near 0 or 1) where Wald intervals fail. Empirical coverage
is 93.4%, consistent with the target range identified by
@brown2001interval.

## Right-Truncation Delay Correction

Sequences collected in recent weeks are incompletely reported because
reporting delays of 1–4 weeks cause right-truncation: sequences with
long delays are systematically missing from recent time points.
[`surv_estimate_delay()`](https://cuiweig.github.io/survinger/reference/surv_estimate_delay.md)
fits parametric delay distributions (negative binomial, Poisson,
lognormal) or a nonparametric EM-based estimator, accounting for
right-truncation in the likelihood.
[`surv_nowcast_lineage()`](https://cuiweig.github.io/survinger/reference/surv_nowcast_lineage.md)
inflates observed counts by the estimated reporting probability to
recover the underlying trajectory.

[`surv_adjusted_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_adjusted_prevalence.md)
combines design weighting and delay correction in a single estimate,
propagating uncertainty from both sources.

## Resource Allocation Optimisation

Given a fixed total sequencing capacity,
[`surv_optimize_allocation()`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md)
solves for the per-stratum allocation minimising one of three
objectives:

- **Minimum MSE**: Neyman allocation, directing more sequencing to
  strata with higher variability.
- **Maximum detection**: Maximising the probability of detecting a
  variant at a specified prevalence.
- **Minimum imbalance**: Minimising deviation from
  population-proportional representation.

A `min_per_stratum` constraint ensures operational feasibility by
preventing any stratum from receiving zero allocation.

## Diagnostic and Reporting Tools

[`surv_report()`](https://cuiweig.github.io/survinger/reference/surv_report.md)
generates a comprehensive system diagnostic including sequencing
inequality (Gini coefficient), effective sample size (Kish’s formula),
detection power, and mean bias of the naive estimator.
[`surv_sensitivity()`](https://cuiweig.github.io/survinger/reference/surv_sensitivity.md)
evaluates robustness by comparing all three estimators on the same data.
[`surv_quality()`](https://cuiweig.github.io/survinger/reference/surv_quality.md)
produces a single-row summary suitable for automated dashboards.

# Example Usage

``` r
library(survinger)
data(sarscov2_surveillance)

# Encode surveillance design
design <- surv_design(
  sarscov2_surveillance,
  strata          = ~region,
  sequencing_rate = sequencing_rate,
  population      = population
)

# Design-weighted prevalence
prev <- surv_lineage_prevalence(design, lineage = "BA.2.86")

# Compare weighted vs naive estimates
surv_compare_estimates(design, lineage = "BA.2.86")

# Fit reporting delay distribution
delay <- surv_estimate_delay(design, distribution = "nbinom")

# Combined design + delay correction
adjusted <- surv_adjusted_prevalence(design, delay, lineage = "BA.2.86")

# Optimise sequencing allocation
alloc <- surv_optimize_allocation(design, objective = "min_mse",
                                  total_capacity = 500,
                                  min_per_stratum = 20)

# System diagnostic
surv_report(design)
```

# Validation

Three levels of validation establish methodological correctness. First,
the Hajek estimator reproduces results from
[`survey::svymean()`](https://rdrr.io/pkg/survey/man/surveysummary.html)
on equivalent designs \[@lumley2004analysis\]. Second, Wilson score
confidence intervals achieve 93.4% empirical coverage. Third, delay
distribution MLE recovers true parameters with 0.5% error at n = 5,000.

The real-data ECDC case study — five countries, 99,093 sequences, Gini
coefficient 0.54 indicating severe sequencing inequality — demonstrates
that design weighting produces a mean 3.8 percentage point correction
relative to naive estimation, with Neyman allocation reducing MSE by 27%
relative to equal allocation.

# Included Data

`sarscov2_surveillance` provides a simulated dataset with five regions,
26 weeks, three sample sources (clinical, wastewater, sentinel), four
lineages, and sequencing rates ranging from 0.5% to 15%. The dataset
includes ground truth prevalences, enabling direct assessment of
estimator bias and variance.

# Design Philosophy

`survinger` follows tidyverse conventions and provides broom-compatible
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) methods.
S3 classes with [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
methods produce publication-quality figures via a custom
[`theme_survinger()`](https://cuiweig.github.io/survinger/reference/theme_survinger.md).
The package uses `checkmate` for input validation, ensuring informative
error messages when surveillance data do not conform to expected
formats.

# Availability

`survinger` is available on CRAN at
<https://CRAN.R-project.org/package=survinger> and on GitHub at
<https://github.com/CuiweiG/survinger>. Four vignettes cover the
quick-start workflow, allocation optimisation, delay correction, and the
ECDC real-data case study.

# References
