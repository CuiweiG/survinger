# survinger and related R tools: scope and boundaries

Last checked: 2026-11-26, for survinger 0.1.1 (the CRAN release of
2026-04-05; the R code has not changed since). Versions and dates of
the other tools were read on that day from CRAN DESCRIPTION files, CRAN
reference manuals, GitHub and r-universe; check them again before
relying on them.

## Summary

survinger estimates pathogen lineage prevalence from sequence-level
records when the fraction of cases that is sequenced differs between
strata such as regions or sample sources. Around one design object it
provides inverse-sequencing-rate weights, per-period prevalence estimates,
a reporting-delay fit with a count nowcast, simple allocation rules, and
detection-probability calculations.

Each of these tasks has a more complete implementation in another package:
'survey' for design-based estimation, several packages for sample
allocation, and several packages for nowcasting. What survinger adds is
the combination for this one data shape, with a small dependency stack
(pure R, no compiled code, no Stan or JAGS). Users who need clustered or
multistage designs, replicate weights, calibration, exact allocation under
bounds, or nowcasts that carry the uncertainty of the delay distribution
should use the specialized tools, alone or next to survinger.

## What survinger implements (0.1.1)

| Component | Functions (file) | What the code does |
|---|---|---|
| Design object | `surv_design()`, `surv_update_rates()`, `surv_set_weights()`, `surv_filter()` (R/01-design.R) | Stores sequence records, the strata, one sequencing rate per stratum and a population table. Weight = 1 / rate (lines 104-109). Weights do not change over time unless the strata include a period column and rates are supplied per stratum and period. |
| Prevalence | `surv_lineage_prevalence()`, `surv_naive_prevalence()` (R/03-prevalence.R); `surv_prevalence_by()`, `surv_sensitivity()` (R/12-advanced.R); `surv_estimate()` (R/13-infrastructure.R) | For each period (ISO week by default): weighted ratio (Hajek) estimate with a stratified with-replacement linearization SE; a "poststratified" option that weights stratum proportions by `pop_total` or `n_positive` shares; Wilson score intervals computed from an effective sample size (Kish for the weighted estimator). The `"horvitz_thompson"` option returns the same point estimate as `"hajek"` (lines 163-182). |
| Delay and nowcast | `surv_estimate_delay()`, `surv_reporting_probability()`, `surv_nowcast_lineage()` (R/04-delay.R) | Negative binomial, Poisson, lognormal or nonparametric delay distribution; the likelihood conditions each delay on the time elapsed since collection, which corrects for right truncation (lines 236-252). The nowcast divides the observed count of a period by the reporting probability at the period's median collection date (lines 190-219). |
| Combined estimate | `surv_adjusted_prevalence()` (R/05-adjusted.R) | Divides stratum-by-period counts by the reporting probability and combines them with the stratum weights; the SE adds a design term and a delay term; normal-approximation interval (lines 59-106). |
| Allocation | `surv_optimize_allocation()`, `surv_compare_allocations()` (R/02-allocation.R) | Closed-form shares with rounding: Neyman-type (`"min_mse"`), proportional (`"min_imbalance"`), and an inverse-prevalence rule (`"max_detection"`, labeled `detection_heuristic` in the output); a minimum per stratum. |
| Detection | `surv_detection_probability()`, `surv_required_sequences()` (R/06-detection.R); `surv_power_curve()` (R/12-advanced.R) | Probability of at least one detection from the average weekly number of sequences (binomial; a Poisson approximation for a threshold above one), assuming the same prevalence in every stratum. |
| Simulation and output | `surv_simulate()` (R/07-simulate.R); `surv_report()` (R/11-report.R); `surv_quality()`, `surv_table()`, `surv_bind()` (R/13-infrastructure.R); plot and tidy methods (R/09-plot-methods.R, R/10-tidy-methods.R) | Synthetic data with unequal rates, lineage dynamics and negative binomial delays; text and table summaries. |

Imports are checkmate, cli, dplyr, generics, ggplot2, methods, purrr,
rlang, stats and tibble. Of the packages in Suggests, nloptr, MASS, Matrix
and vctrs are not called from R/, and 'survey' is used only in tests and
documentation.

## Comparable tools

Release information is the latest version on CRAN unless stated
otherwise. Function-level statements come from each package's reference
manual.

### (i) Design-based prevalence estimation

| Tool | Relevant scope | Latest release |
|---|---|---|
| 'survey' (Thomas Lumley) | Stratified, clustered, multistage and unequal-probability designs, with or without finite population corrections; variances by Taylor linearization or replicate weights (`as.svrepdesign()`: jackknife, BRR, bootstrap, Fay); post-stratification (`postStratify()`), raking (`rake()`) and calibration (`calibrate()`); domain estimates (`svyby()`, `subset()`); intervals for proportions (`svyciprop()`: logit, likelihood, asin, beta (Korn and Graubard), mean, xlogit, wilson). Lumley (2004). | 4.5, 2026-02-24 ([CRAN](https://CRAN.R-project.org/package=survey), [manual](https://cran.r-project.org/web/packages/survey/refman/survey.html)) |
| 'srvyr' | dplyr-style interface that imports 'survey' and wraps its functions; `survey_mean(proportion = TRUE)` calls `svyciprop()`. | 1.3.1, 2026-02-02 ([CRAN](https://CRAN.R-project.org/package=srvyr), [manual](https://cran.r-project.org/web/packages/srvyr/refman/srvyr.html)) |
| 'sampling' (Yves Tillé, Alina Matei) | Sample selection (stratified, cluster, multistage, unequal-probability and balanced designs) and estimators: `HTestimator()`, `Hajekestimator()`, calibration (`calib()`), post-stratification (`postest()`, `poststrata()`), variance estimators (`varHT()`, `varest()`, `vartaylor_ratio()`). I found no allocation function in the 2.11 manual. | 2.11, 2025-07-10 ([CRAN](https://CRAN.R-project.org/package=sampling), [manual](https://cran.r-project.org/web/packages/sampling/refman/sampling.html)) |

### (ii) Sample allocation

| Tool | Relevant scope | Latest release |
|---|---|---|
| 'PracTools' | `strAlloc()`: proportional, Neyman, cost-constrained and variance-constrained allocation for stratified simple random samples; sample sizes for one- to three-stage designs. Companion to Valliant, Dever and Kreuter (2018). | 1.8, 2026-11-05 ([CRAN](https://CRAN.R-project.org/package=PracTools), [manual](https://cran.r-project.org/web/packages/PracTools/refman/PracTools.html)) |
| 'optimall' | `optimum_allocation()`: Neyman or Wright exact integer allocation, with lower and upper bounds under the WrightIII method; multiwave and multiphase designs. Yang et al. (2025). | 1.4.0, 2026-04-03 ([CRAN](https://CRAN.R-project.org/package=optimall), [manual](https://cran.r-project.org/web/packages/optimall/refman/optimall.html)) |
| 'stratallo' | Exact algorithms: Neyman-Tschuprow allocation with lower and upper bounds (`opt()`), minimum cost under a variance constraint (`optcost()`), multidomain allocation (`dopt()`), rounding helpers. | 3.0.1, 2026-03-12 ([CRAN](https://CRAN.R-project.org/package=stratallo), [manual](https://cran.r-project.org/web/packages/stratallo/refman/stratallo.html)) |
| 'MOSAlloc' | Multiobjective allocation with precision and cost constraints at the subpopulation level (`mosalloc()`, `mosallocSTRS()`), single- and multistage designs. | 1.2.6, 2026-06-07 ([CRAN](https://CRAN.R-project.org/package=MOSAlloc), [manual](https://cran.r-project.org/web/packages/MOSAlloc/refman/MOSAlloc.html)) |
| 'SamplingStrata' | Optimal stratification of a sampling frame, with Bethel-Chromy allocation to meet multivariate precision targets. | 1.5-5, 2025-10-08 ([CRAN](https://CRAN.R-project.org/package=SamplingStrata)) |
| 'allocation' | Exact optimal allocation methods of Wright (2012, 2017). | 0.1.0, 2025-08-27 ([CRAN](https://CRAN.R-project.org/package=allocation)) |

### (iii) Reporting-delay nowcasting

| Tool | Relevant scope | Latest release |
|---|---|---|
| `surveillance::nowcast()` | Takes one row per case with event and report dates, the same input shape survinger uses. Methods include a Bayesian method under right truncation (`"bayes.trunc"`), the method of Lawless (1994) and a variant with change points in the delay distribution (`"bayes.trunc.ddcp"`). Höhle and an der Heiden (2014). | 'surveillance' 1.26.1, 2026-07-30 ([CRAN](https://CRAN.R-project.org/package=surveillance), [manual](https://cran.r-project.org/web/packages/surveillance/refman/surveillance.html)) |
| 'NobBS' | Bayesian nowcast that estimates the delay distribution and the epidemic curve together; stratified version `NobBS.strat()`; requires JAGS. McGough et al. (2020). | 1.1.1, 2026-06-01 ([CRAN](https://CRAN.R-project.org/package=NobBS), [manual](https://cran.r-project.org/web/packages/NobBS/refman/NobBS.html)) |
| 'epinowcast' | Modular Bayesian framework for nowcasting and delay estimation: parametric or nonparametric delays, time-varying effects, partial pooling across strata, missing data. Model fitting needs 'cmdstanr' and CmdStan, which are not installed automatically ([README](https://github.com/epinowcast/epinowcast/blob/v0.7.0/README.md)). Not on CRAN. | v0.7.0, GitHub release 2026-07-11 ([release](https://github.com/epinowcast/epinowcast/releases/tag/v0.7.0), [r-universe](https://epinowcast.r-universe.dev/api/packages/epinowcast)) |
| 'EpiNow2' | Reproduction-number estimation with a renewal model in Stan; `estimate_truncation()` estimates a truncation distribution from successive data snapshots and passes its uncertainty to `epinow()` and `estimate_infections()`. | 1.9.0, 2026-06-16 ([CRAN](https://CRAN.R-project.org/package=EpiNow2), [manual](https://cran.r-project.org/web/packages/EpiNow2/refman/EpiNow2.html)) |
| 'baselinenowcast' | Nowcasts from empirical delay distributions with uncertainty estimated from past performance; no Stan; intended as a baseline method. Johnson et al. (Wellcome Open Research). | 0.2.2, 2026-11-04 ([CRAN](https://CRAN.R-project.org/package=baselinenowcast), [manual](https://cran.r-project.org/web/packages/baselinenowcast/refman/baselinenowcast.html)) |
| 'nowcastr' | Chain-ladder nowcasting, with delay-based or completeness-curve estimation. | 0.2.1, 2026-07-07 ([CRAN](https://CRAN.R-project.org/package=nowcastr)) |

### (iv) Genomic-surveillance sampling and variant proportions

| Tool | Relevant scope | Latest release |
|---|---|---|
| 'phylosamp' | Sample size and detection probability for variant surveillance (`vartrack_samplesize_detect()`, `vartrack_samplesize_prev()`, `vartrack_prob_detect()`), with a coefficient-of-detection ratio for differences in detection between variants. Wohl et al. (2023). | 1.0.1, 2023-05-23 ([CRAN](https://CRAN.R-project.org/package=phylosamp), [manual](https://cran.r-project.org/web/packages/phylosamp/refman/phylosamp.html)); last push to [GitHub](https://github.com/HopkinsIDD/phylosamp) 2023-05-23 |
| CDC SARS-CoV-2 genomic surveillance code | R scripts, not a package. Per the repository README, the methods produce the variant proportions on CDC's COVID Data Tracker. Survey design with laboratory source as clusters within state-by-week strata and weights from the estimated number of infections each sequence represents; proportions from `survey::svyciprop()` (xlogit, or a modified Korn and Graubard interval); survey-weighted multinomial regression (`nnet::multinom()` with variance from `survey::svyrecvar()`) for the most recent periods. Lambrou et al. (2022). | No releases; last commit on the default branch 2025-07-09 ([repository](https://github.com/CDCgov/SARS-CoV-2_Genomic_Surveillance)) |
| 'lineagefreq' (same maintainer as survinger) | Lineage frequency dynamics from counts: multinomial logistic regression, hierarchical partial pooling, Piantham approximation for relative reproduction numbers, backtesting, `sequencing_power()`. | 0.2.0, 2026-04-03 ([CRAN](https://CRAN.R-project.org/package=lineagefreq), [manual](https://cran.r-project.org/web/packages/lineagefreq/refman/lineagefreq.html)) |

## Boundaries

### What survinger adds

- One object for sequence-level records (collection date, report date,
  lineage, stratum, stratum sequencing rate, population table) that the
  weighting, delay and allocation functions all accept (`surv_design()`,
  R/01-design.R). 'survey' expects the weights or selection probabilities
  as a column; survinger derives them from stratum sequencing rates given
  as a formula, a named vector or a data frame.
- Per-period lineage prevalence next to an unweighted estimate, a
  design-effect table and a method comparison in a few calls
  (`surv_lineage_prevalence()`, `surv_naive_prevalence()`,
  `surv_design_effect()`, `surv_sensitivity()`).
- Planning helpers that work on the same object: allocation rules and
  detection probabilities (`surv_optimize_allocation()`,
  `surv_detection_probability()`, `surv_required_sequences()`,
  `surv_power_curve()`).
- A simulator with known truth for unequal sequencing rates and negative
  binomial reporting delays (`surv_simulate()`, R/07-simulate.R), useful
  for teaching and for testing pipelines.
- A small dependency stack: pure R, no compiled code, no Stan or JAGS.

The estimators are standard: the Hajek ratio with its linearization
variance, Neyman allocation, the Wilson score interval, and division of
counts by a reporting probability. The package does not add new
estimators.

### What other tools do better or more completely

- Designs. 'survey' handles clusters (for example laboratories as primary
  sampling units, as in the CDC design), multistage and PPS designs, and
  finite population corrections. survinger handles one-stage
  stratification with known rates and no finite population correction in
  the Hajek SE.
- Variance estimation. 'survey' offers linearization, replicate weights
  and design degrees of freedom. survinger offers one linearization
  formula for the Hajek estimator.
- Intervals. `svyciprop()` offers several interval methods built on the
  design-based variance. survinger's Wilson intervals use the Kish
  effective sample size, which reflects the spread of the weights but not
  the effect of stratification on the variance; the reported `se` column
  is not used for the interval.
- Weight adjustment. 'survey' post-stratifies, rakes and calibrates to
  known totals (for example case counts by age group); 'sampling' provides
  calibration weights and post-stratified estimators. survinger takes the
  rates as given; its "poststratified" option weights the design strata
  by population shares and does not adjust to other margins.
- Domains. `svyby()` and `subset()` in 'survey' keep the full design for
  domain variances. `surv_prevalence_by()` refits a filtered design
  (R/12-advanced.R:45-60), so domain sample sizes are treated as fixed.
- Several lineages at once. Multinomial models (CDC code, 'lineagefreq')
  estimate all lineage proportions jointly and smooth over time. survinger
  estimates one lineage at a time, period by period.
- Allocation. 'optimall', 'stratallo', 'MOSAlloc', 'PracTools' and
  'SamplingStrata' solve allocation problems exactly or with explicit
  upper bounds, costs, precision targets or several domains. survinger
  rounds closed-form shares and enforces only a minimum per stratum.
- Nowcasting. The Bayesian tools ('surveillance' with `"bayes.trunc"`,
  'NobBS', 'epinowcast', 'EpiNow2') carry the uncertainty of the
  estimated delay distribution into the nowcast, and several of them
  allow delays or incidence to change over time; 'baselinenowcast'
  derives nowcast uncertainty from past nowcast errors. survinger divides
  counts by a point estimate of the reporting probability.
- Sample size. 'phylosamp' allows differences in detection between
  variants. survinger assumes equal detection and equal prevalence across
  strata.

## When to use which

Use 'survey' (or 'srvyr') instead of survinger when:

- sequences are clustered (laboratories, sites, outbreaks) or selected in
  more than one stage;
- variances from replicate weights, finite population corrections or
  design degrees of freedom are needed;
- weights must be calibrated or raked to known totals, for example case
  counts by age, sex and region;
- estimates for domains that cut across strata are needed;
- results feed official statistics and must follow an established survey
  estimation standard.

Use 'survey' together with survinger when the survinger design object is
convenient for data handling but design-based SEs and intervals are
needed, or to cross-check survinger's estimates. This sketch mirrors
the package's own test in tests/testthat/test-03-prevalence.R:71-85;
it runs with survinger 0.1.1 and 'survey' 4.5:

```r
library(survinger)
sim <- surv_simulate(n_regions = 5, n_weeks = 12, seed = 1)
d <- surv_design(sim$sequences, ~ region,
                 sim$population[c("region", "seq_rate")], sim$population)
dat <- merge(d$data, d$strata_info[c("region", "seq_rate")], by = "region")
dat$wt <- 1 / dat$seq_rate
des <- survey::svydesign(ids = ~1, strata = ~region, weights = ~wt,
                         data = dat)
survey::svyciprop(~ I(lineage == "BA.2.86"),
                  subset(des, epiweek == "2024-W10"), method = "beta")
```

For other tasks:

- Allocation: use 'optimall', 'stratallo', 'MOSAlloc' or 'PracTools'
  when the allocation must respect upper bounds (for example the number of
  samples available per stratum), costs, precision targets or several
  domains. survinger's rules give a quick first split.
- Nowcasting counts: use `surveillance::nowcast()`, 'NobBS',
  'epinowcast' or 'baselinenowcast' when the uncertainty of the delay
  distribution or changes in delays over time matter.
- Variant proportions in the most recent periods: use a multinomial model
  (the CDC code; 'lineagefreq'), which borrows strength across time and
  lineages.
- Sample size with differences in detection between variants: use
  'phylosamp'.

survinger is a reasonable choice when sequencing is stratified in one
stage with known sequencing fractions per stratum, the data are
sequence-level records, and the goal is a quick weekly weighted estimate
with an unweighted comparison, simple allocation and detection planning,
or simulated test data, without heavy dependencies.

## Sources

Accessed 2026-11-26.

CRAN package pages (version and publication date from each package's
DESCRIPTION file):
[survey](https://CRAN.R-project.org/package=survey),
[srvyr](https://CRAN.R-project.org/package=srvyr),
[sampling](https://CRAN.R-project.org/package=sampling),
[PracTools](https://CRAN.R-project.org/package=PracTools),
[optimall](https://CRAN.R-project.org/package=optimall),
[stratallo](https://CRAN.R-project.org/package=stratallo),
[MOSAlloc](https://CRAN.R-project.org/package=MOSAlloc),
[SamplingStrata](https://CRAN.R-project.org/package=SamplingStrata),
[allocation](https://CRAN.R-project.org/package=allocation),
[surveillance](https://CRAN.R-project.org/package=surveillance),
[NobBS](https://CRAN.R-project.org/package=NobBS),
[EpiNow2](https://CRAN.R-project.org/package=EpiNow2),
[baselinenowcast](https://CRAN.R-project.org/package=baselinenowcast),
[nowcastr](https://CRAN.R-project.org/package=nowcastr),
[phylosamp](https://CRAN.R-project.org/package=phylosamp),
[lineagefreq](https://CRAN.R-project.org/package=lineagefreq),
[survinger](https://CRAN.R-project.org/package=survinger).

epinowcast is not on CRAN: no page at
<https://cran.r-project.org/web/packages/epinowcast/> and no CRAN
archive. Release:
<https://github.com/epinowcast/epinowcast/releases/tag/v0.7.0>;
distribution: <https://epinowcast.r-universe.dev/api/packages/epinowcast>.

CDC code: <https://github.com/CDCgov/SARS-CoV-2_Genomic_Surveillance>
(README and `Code/variant surveillance system.R`); commit history from
<https://api.github.com/repos/CDCgov/SARS-CoV-2_Genomic_Surveillance/commits>.

Papers (DOIs checked with Crossref):

- Lumley T (2004). Analysis of complex survey samples. Journal of
  Statistical Software 9(8). <https://doi.org/10.18637/jss.v009.i08>
- Höhle M, an der Heiden M (2014). Bayesian nowcasting during the STEC
  O104:H4 outbreak in Germany, 2011. Biometrics 70(4):993-1002.
  <https://doi.org/10.1111/biom.12194>
- McGough SF, Johansson MA, Lipsitch M, Menzies NA (2020). Nowcasting by
  Bayesian Smoothing. PLOS Computational Biology 16(4):e1007735.
  <https://doi.org/10.1371/journal.pcbi.1007735>
- Yang JB, Lumley T, Shepherd BE, Shaw PA (2025). Optimum allocation for
  adaptive multi-wave sampling in R: the R package optimall. Journal of
  Statistical Software 114(10). <https://doi.org/10.18637/jss.v114.i10>
- Wohl S, Lee EC, DiPrete BL, Lessler J (2023). Sample size calculations
  for pathogen variant surveillance in the presence of biological and
  systematic biases. Cell Reports Medicine 4(5):101022.
  <https://doi.org/10.1016/j.xcrm.2023.101022>
- Lambrou AS, Shirk P, Steele MK, et al. (2022). Genomic surveillance for
  SARS-CoV-2 variants: predominance of the Delta (B.1.617.2) and Omicron
  (B.1.1.529) variants, United States, June 2021 to January 2022. MMWR
  71(6):206-211. <https://doi.org/10.15585/mmwr.mm7106a4>
- Johnson K, et al. Baseline nowcasting methods for handling delays in
  epidemiological data. Wellcome Open Research 10:614, version 2.
  <https://doi.org/10.12688/wellcomeopenres.25027.2>

## Not verified

- Agreement of survinger's Hajek SE with `survey::svymean()`: the SE
  equals the standard with-replacement linearization formula in a manual
  check; the outputs of the two packages were not compared.
- Whether the CDC code in its current form is what produces the COVID Data
  Tracker estimates: taken from the repository README.
- Maintenance beyond the latest release date (open issues, CRAN check
  results, response times) was not assessed for any tool.
