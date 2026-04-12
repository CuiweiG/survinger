# When Counting Sequences Is Not Enough: Survey Sampling Methods for Genomic Surveillance

Denmark sequences 12% of its positive SARS-CoV-2 cases. Romania
sequences 0.3%. If you estimate lineage prevalence across Europe by
simply counting sequences, you will get Denmark’s epidemic, not
Europe’s. On real European Centre for Disease Prevention and Control
(ECDC) data encompassing 99,093 sequences from five countries, this
naive approach deviates from design-corrected estimates by up to 14
percentage points.

This is not a novel statistical insight. Survey sampling theory solved
the problem of unequal inclusion probabilities decades ago — Horvitz and
Thompson published their foundational estimator in 1952. Yet the genomic
surveillance community has been slow to adopt these methods, in part
because no off-the-shelf tool translated survey sampling concepts into
the language and data structures of pathogen genomics.

The `survinger` R package bridges that gap. It implements
design-adjusted prevalence estimation, right-truncation-corrected delay
distributions, and Neyman-optimal resource allocation within a workflow
designed specifically for genomic surveillance operations. It is now
available on CRAN.

## The Three Problems

Genomic surveillance data suffer from three distinct sources of bias,
each requiring its own correction.

**Unequal sequencing rates.** Laboratories, regions, and countries
sequence at vastly different rates. Without inverse-probability
weighting, prevalence estimates are dominated by high-sequencing
jurisdictions regardless of their population size. The Hajek estimator —
a ratio form of inverse-probability weighting — corrects this while
maintaining efficiency superior to the Horvitz-Thompson estimator.

**Reporting delays.** Sequences collected today may not appear in public
databases for one to four weeks. This right-truncation means that recent
time points always appear artificially low, creating a systematic
downward bias in the most operationally relevant estimates. Maximum
likelihood fitting of truncated delay distributions (negative binomial,
Poisson, lognormal, or nonparametric) recovers the underlying delay
structure and enables nowcasting.

**Inefficient allocation.** Fixed sequencing budgets must be distributed
across strata. Equal allocation ignores heterogeneity in population size
and variant prevalence. Proportional allocation ignores differences in
variability. Neyman allocation minimises estimation variance by
directing more sequencing effort toward strata with higher variability —
achieving, in the ECDC case study, a 27% reduction in mean squared error
relative to equal allocation.

## A Worked Example

The workflow begins by encoding the surveillance design — strata,
sequencing rates, and population denominators:

``` r
library(survinger)
data(sarscov2_surveillance)

design <- surv_design(
  sarscov2_surveillance,
  strata          = ~region,
  sequencing_rate = sequencing_rate,
  population      = population
)
```

Design-weighted prevalence estimation for a target lineage uses the
Hajek estimator by default, with Wilson score confidence intervals that
maintain proper coverage even at extreme proportions:

``` r
prev <- surv_lineage_prevalence(design, lineage = "BA.2.86")
plot(prev)
```

To see how much the weighting matters,
[`surv_compare_estimates()`](https://cuiweig.github.io/survinger/reference/surv_compare_estimates.md)
places the naive and design-adjusted estimates side by side:

``` r
surv_compare_estimates(design, lineage = "BA.2.86")
```

On the included simulation dataset — five regions with sequencing rates
ranging from 0.5% to 15% — the gap between weighted and unweighted
estimates is immediately apparent and varies by week as the lineage
dynamics interact with the sequencing inequality.

## Delay Correction and Nowcasting

Fitting the delay distribution accounts for right-truncation: sequences
with long delays are systematically missing from recent time points, so
the observed delay distribution is biased toward short delays. The
likelihood is modified accordingly:

``` r
delay <- surv_estimate_delay(design, distribution = "nbinom")
surv_reporting_probability(delay, days_ago = 7)
```

Nowcasting inflates observed counts by the reporting probability to
recover the true underlying trajectory:

``` r
nowcast <- surv_nowcast_lineage(design, delay, lineage = "BA.2.86")
plot(nowcast)
```

The combined correction — design weighting plus delay adjustment — is
available in a single function that propagates uncertainty from both
sources:

``` r
adjusted <- surv_adjusted_prevalence(design, delay, lineage = "BA.2.86")
```

## Optimising Sequencing Resources

Given a fixed total sequencing capacity,
[`surv_optimize_allocation()`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md)
solves for the per-stratum allocation that minimises a chosen objective:

``` r
alloc <- surv_optimize_allocation(
  design,
  objective       = "min_mse",
  total_capacity  = 500,
  min_per_stratum = 20
)
plot(alloc)
```

Three objectives are available: minimise prevalence MSE (Neyman
allocation), maximise detection probability for a rare variant, or
minimise representational imbalance across strata. The `min_per_stratum`
floor ensures that no region is entirely excluded, a practical
constraint that pure optimisation would otherwise violate.

[`surv_compare_allocations()`](https://cuiweig.github.io/survinger/reference/surv_compare_allocations.md)
benchmarks all strategies against equal and proportional baselines,
providing the evidence base for resource allocation decisions.

## System Diagnostics

[`surv_report()`](https://cuiweig.github.io/survinger/reference/surv_report.md)
generates a comprehensive diagnostic summary: total sequences, number of
strata, date range, sequencing inequality (Gini coefficient), effective
sample size (Kish’s formula), detection power for a specified
prevalence, and mean bias of the naive estimator.
[`surv_quality()`](https://cuiweig.github.io/survinger/reference/surv_quality.md)
distils these into a single-row tibble suitable for automated monitoring
dashboards.

For methodological robustness,
[`surv_sensitivity()`](https://cuiweig.github.io/survinger/reference/surv_sensitivity.md)
runs all three estimators (Hajek, Horvitz-Thompson, post-stratified)
plus a delay-adjusted variant on the same data and reports how sensitive
the conclusions are to estimator choice.

## Validation

The package has been validated at three levels. First, the Hajek
estimator produces results that match
[`survey::svymean()`](https://rdrr.io/pkg/survey/man/surveysummary.html)
(Lumley 2004) on equivalent designs. Second, Wilson score intervals
achieve 93.4% empirical coverage, consistent with the target range
identified by Brown, Cai, and DasGupta (2001). Third, the delay
distribution MLE recovers true parameters with 0.5% error at sample
sizes of 5,000 — adequate for most national surveillance systems.

The real-data ECDC case study (99,093 sequences, five countries, Gini
coefficient 0.54) demonstrates that design weighting produces a mean 3.8
percentage point correction relative to naive estimation — a difference
that is operationally meaningful for variant risk assessments and
resource allocation decisions.

## Intended Audience

`survinger` is aimed at epidemiologists, genomic surveillance
coordinators, and public health statisticians who work with
variant-resolved sequence data. It does not require familiarity with
formal survey sampling theory; the functions abstract the weighting
mechanics behind an interface expressed in surveillance-native
terminology (lineages, sequencing rates, reporting delays).

The package follows tidyverse conventions, provides broom-compatible
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`glance()`](https://generics.r-lib.org/reference/glance.html) methods,
and includes publication-quality plotting via a custom
[`theme_survinger()`](https://cuiweig.github.io/survinger/reference/theme_survinger.md).
Four vignettes cover the quick-start workflow, allocation optimisation,
delay correction, and the ECDC real-data case study.

``` r
install.packages("survinger")
```

Source code and documentation are on
[GitHub](https://github.com/CuiweiG/survinger).

Unequal sequencing is not going away. As genomic surveillance matures
from pandemic response to routine public health infrastructure, the
tools applied to its data must reflect the sampling design under which
that data was collected. Classical survey methods, adapted to the
genomic context, provide a rigorous and well-understood foundation for
that work.

------------------------------------------------------------------------

*Cuiwei Gao is a health data analyst and R developer. `survinger` is
available on [CRAN](https://CRAN.R-project.org/package=survinger).*
