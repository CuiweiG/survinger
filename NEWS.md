# survinger 0.1.1

* Fix test failure on M1 Mac (ARM64): Wilson CI is not centred
  on the point estimate by construction. Added floating-point
  tolerance to the CI-containment test (test-03-prevalence.R:63).
  Only the test was affected; no changes to estimation code.

# survinger 0.1.0

* Initial release.
* Core functions: `surv_design()`, `surv_optimize_allocation()`,
  `surv_lineage_prevalence()`, `surv_estimate_delay()`,
  `surv_nowcast_lineage()`, `surv_adjusted_prevalence()`,
  `surv_detection_probability()`, `surv_required_sequences()`.
* Three prevalence estimators: Horvitz-Thompson, Hajek, post-stratified.
* Three allocation objectives: min_mse, max_detection, min_imbalance.
* Reporting delay estimation with right-truncation correction.
* Diagnostic tools: `surv_compare_estimates()`, `surv_design_effect()`.
* Example dataset: `sarscov2_surveillance`.
