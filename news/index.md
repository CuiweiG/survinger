# Changelog

## survinger 0.1.1

CRAN release: 2026-04-05

- Fix test failure on M1 Mac (ARM64): Wilson CI is not centred on the
  point estimate by construction. Added floating-point tolerance to the
  CI-containment test (test-03-prevalence.R:63). Only the test was
  affected; no changes to estimation code.

## survinger 0.1.0

CRAN release: 2026-04-02

- Initial release.
- Core functions:
  [`surv_design()`](https://cuiweig.github.io/survinger/reference/surv_design.md),
  [`surv_optimize_allocation()`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md),
  [`surv_lineage_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md),
  [`surv_estimate_delay()`](https://cuiweig.github.io/survinger/reference/surv_estimate_delay.md),
  [`surv_nowcast_lineage()`](https://cuiweig.github.io/survinger/reference/surv_nowcast_lineage.md),
  [`surv_adjusted_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_adjusted_prevalence.md),
  [`surv_detection_probability()`](https://cuiweig.github.io/survinger/reference/surv_detection_probability.md),
  [`surv_required_sequences()`](https://cuiweig.github.io/survinger/reference/surv_required_sequences.md).
- Three prevalence estimators: Horvitz-Thompson, Hajek, post-stratified.
- Three allocation objectives: min_mse, max_detection, min_imbalance.
- Reporting delay estimation with right-truncation correction.
- Diagnostic tools:
  [`surv_compare_estimates()`](https://cuiweig.github.io/survinger/reference/surv_compare_estimates.md),
  [`surv_design_effect()`](https://cuiweig.github.io/survinger/reference/surv_design_effect.md).
- Example dataset: `sarscov2_surveillance`.
