# Package index

## Design construction

Create and modify surveillance designs with inverse-probability weights.

- [`print(`*`<surv_design>`*`)`](https://cuiweig.github.io/survinger/reference/surv_design.md)
  [`summary(`*`<surv_design>`*`)`](https://cuiweig.github.io/survinger/reference/surv_design.md)
  [`print(`*`<summary.surv_design>`*`)`](https://cuiweig.github.io/survinger/reference/surv_design.md)
  [`surv_design()`](https://cuiweig.github.io/survinger/reference/surv_design.md)
  : Create a genomic surveillance design object
- [`surv_update_rates()`](https://cuiweig.github.io/survinger/reference/surv_update_rates.md)
  : Update sequencing rates in a surveillance design
- [`surv_set_weights()`](https://cuiweig.github.io/survinger/reference/surv_set_weights.md)
  : Override design weights with custom values
- [`surv_filter()`](https://cuiweig.github.io/survinger/reference/surv_filter.md)
  : Subset a surveillance design by filter criteria

## Prevalence estimation

Weighted and unweighted lineage prevalence estimators.

- [`print(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md)
  [`as.data.frame(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md)
  [`surv_lineage_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md)
  : Estimate lineage prevalence with design weights
- [`surv_naive_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_naive_prevalence.md)
  : Compute naive (unweighted) lineage prevalence
- [`surv_prevalence_by()`](https://cuiweig.github.io/survinger/reference/surv_prevalence_by.md)
  : Estimate prevalence by subgroup
- [`surv_estimate()`](https://cuiweig.github.io/survinger/reference/surv_estimate.md)
  : Pipe-friendly surveillance analysis

## Delay correction & nowcasting

Right-truncation-corrected delay fitting and delay-adjusted nowcasts.

- [`print(`*`<surv_delay_fit>`*`)`](https://cuiweig.github.io/survinger/reference/surv_estimate_delay.md)
  [`surv_estimate_delay()`](https://cuiweig.github.io/survinger/reference/surv_estimate_delay.md)
  : Estimate reporting delay distribution
- [`surv_reporting_probability()`](https://cuiweig.github.io/survinger/reference/surv_reporting_probability.md)
  : Compute cumulative reporting probability
- [`print(`*`<surv_nowcast>`*`)`](https://cuiweig.github.io/survinger/reference/surv_nowcast_lineage.md)
  [`as.data.frame(`*`<surv_nowcast>`*`)`](https://cuiweig.github.io/survinger/reference/surv_nowcast_lineage.md)
  [`surv_nowcast_lineage()`](https://cuiweig.github.io/survinger/reference/surv_nowcast_lineage.md)
  : Nowcast lineage counts correcting for reporting delays
- [`print(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/surv_adjusted_prevalence.md)
  [`as.data.frame(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/surv_adjusted_prevalence.md)
  [`surv_adjusted_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_adjusted_prevalence.md)
  : Combined design-weighted and delay-adjusted prevalence

## Resource allocation

Optimise sequencing budgets across strata.

- [`print(`*`<surv_allocation>`*`)`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md)
  [`as.data.frame(`*`<surv_allocation>`*`)`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md)
  [`surv_optimize_allocation()`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md)
  : Optimize sequencing allocation across strata
- [`surv_compare_allocations()`](https://cuiweig.github.io/survinger/reference/surv_compare_allocations.md)
  : Compare multiple allocation strategies
- [`surv_required_sequences()`](https://cuiweig.github.io/survinger/reference/surv_required_sequences.md)
  : Required sequences for target detection probability

## Diagnostics & reporting

Power analysis, design effects, sensitivity, and one-page reports.

- [`surv_detection_probability()`](https://cuiweig.github.io/survinger/reference/surv_detection_probability.md)
  : Variant detection probability under current design
- [`surv_power_curve()`](https://cuiweig.github.io/survinger/reference/surv_power_curve.md)
  [`plot(`*`<surv_power_curve>`*`)`](https://cuiweig.github.io/survinger/reference/surv_power_curve.md)
  : Compute power curve for detection across prevalence range
- [`surv_design_effect()`](https://cuiweig.github.io/survinger/reference/surv_design_effect.md)
  : Compute design effect over time
- [`surv_sensitivity()`](https://cuiweig.github.io/survinger/reference/surv_sensitivity.md)
  : Sensitivity analysis across methods
- [`surv_report()`](https://cuiweig.github.io/survinger/reference/surv_report.md)
  : Generate a comprehensive surveillance system report
- [`surv_quality()`](https://cuiweig.github.io/survinger/reference/surv_quality.md)
  : Compute surveillance quality metrics

## Visualisation

Plot methods and comparison helpers.

- [`surv_compare_estimates()`](https://cuiweig.github.io/survinger/reference/surv_compare_estimates.md)
  : Compare weighted vs naive prevalence estimates
- [`surv_plot_sequencing_rates()`](https://cuiweig.github.io/survinger/reference/surv_plot_sequencing_rates.md)
  : Plot sequencing rate inequality across strata
- [`surv_plot_allocation()`](https://cuiweig.github.io/survinger/reference/surv_plot_allocation.md)
  : Plot allocation plan
- [`theme_survinger()`](https://cuiweig.github.io/survinger/reference/theme_survinger.md)
  : Publication-quality ggplot2 theme
- [`plot(`*`<surv_design>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_allocation>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_delay_fit>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_nowcast>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  : Plot methods for survinger objects
- [`surv_power_curve()`](https://cuiweig.github.io/survinger/reference/surv_power_curve.md)
  [`plot(`*`<surv_power_curve>`*`)`](https://cuiweig.github.io/survinger/reference/surv_power_curve.md)
  : Compute power curve for detection across prevalence range

## Tidyverse integration

Broom-style tidying and table helpers.

- [`tidy(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  [`tidy(`*`<surv_nowcast>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  [`tidy(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  [`tidy(`*`<surv_allocation>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  [`tidy(`*`<surv_delay_fit>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  : Extract tidy estimates from survinger objects
- [`glance(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/glance.surv.md)
  [`glance(`*`<surv_delay_fit>`*`)`](https://cuiweig.github.io/survinger/reference/glance.surv.md)
  [`glance(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/glance.surv.md)
  : One-row summary of survinger model
- [`surv_bind()`](https://cuiweig.github.io/survinger/reference/surv_bind.md)
  : Combine multiple prevalence estimates
- [`surv_table()`](https://cuiweig.github.io/survinger/reference/surv_table.md)
  : Format prevalence results for knitr tables

## Simulation

Generate synthetic surveillance datasets for testing and benchmarking.

- [`surv_simulate()`](https://cuiweig.github.io/survinger/reference/surv_simulate.md)
  : Simulate genomic surveillance data

## Package

Package-level documentation.

- [`survinger`](https://cuiweig.github.io/survinger/reference/survinger-package.md)
  [`survinger-package`](https://cuiweig.github.io/survinger/reference/survinger-package.md)
  : survinger: Design-Adjusted Inference for Pathogen Lineage
  Surveillance

## Other

- [`glance(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/glance.surv.md)
  [`glance(`*`<surv_delay_fit>`*`)`](https://cuiweig.github.io/survinger/reference/glance.surv.md)
  [`glance(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/glance.surv.md)
  : One-row summary of survinger model
- [`plot(`*`<surv_design>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_allocation>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_delay_fit>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_nowcast>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  : Plot methods for survinger objects
- [`sarscov2_surveillance`](https://cuiweig.github.io/survinger/reference/sarscov2_surveillance.md)
  : Example SARS-CoV-2 genomic surveillance data
- [`print(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/surv_adjusted_prevalence.md)
  [`as.data.frame(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/surv_adjusted_prevalence.md)
  [`surv_adjusted_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_adjusted_prevalence.md)
  : Combined design-weighted and delay-adjusted prevalence
- [`surv_bind()`](https://cuiweig.github.io/survinger/reference/surv_bind.md)
  : Combine multiple prevalence estimates
- [`surv_compare_allocations()`](https://cuiweig.github.io/survinger/reference/surv_compare_allocations.md)
  : Compare multiple allocation strategies
- [`surv_compare_estimates()`](https://cuiweig.github.io/survinger/reference/surv_compare_estimates.md)
  : Compare weighted vs naive prevalence estimates
- [`print(`*`<surv_design>`*`)`](https://cuiweig.github.io/survinger/reference/surv_design.md)
  [`summary(`*`<surv_design>`*`)`](https://cuiweig.github.io/survinger/reference/surv_design.md)
  [`print(`*`<summary.surv_design>`*`)`](https://cuiweig.github.io/survinger/reference/surv_design.md)
  [`surv_design()`](https://cuiweig.github.io/survinger/reference/surv_design.md)
  : Create a genomic surveillance design object
- [`surv_design_effect()`](https://cuiweig.github.io/survinger/reference/surv_design_effect.md)
  : Compute design effect over time
- [`surv_detection_probability()`](https://cuiweig.github.io/survinger/reference/surv_detection_probability.md)
  : Variant detection probability under current design
- [`surv_estimate()`](https://cuiweig.github.io/survinger/reference/surv_estimate.md)
  : Pipe-friendly surveillance analysis
- [`print(`*`<surv_delay_fit>`*`)`](https://cuiweig.github.io/survinger/reference/surv_estimate_delay.md)
  [`surv_estimate_delay()`](https://cuiweig.github.io/survinger/reference/surv_estimate_delay.md)
  : Estimate reporting delay distribution
- [`surv_filter()`](https://cuiweig.github.io/survinger/reference/surv_filter.md)
  : Subset a surveillance design by filter criteria
- [`print(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md)
  [`as.data.frame(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md)
  [`surv_lineage_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md)
  : Estimate lineage prevalence with design weights
- [`surv_naive_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_naive_prevalence.md)
  : Compute naive (unweighted) lineage prevalence
- [`print(`*`<surv_nowcast>`*`)`](https://cuiweig.github.io/survinger/reference/surv_nowcast_lineage.md)
  [`as.data.frame(`*`<surv_nowcast>`*`)`](https://cuiweig.github.io/survinger/reference/surv_nowcast_lineage.md)
  [`surv_nowcast_lineage()`](https://cuiweig.github.io/survinger/reference/surv_nowcast_lineage.md)
  : Nowcast lineage counts correcting for reporting delays
- [`print(`*`<surv_allocation>`*`)`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md)
  [`as.data.frame(`*`<surv_allocation>`*`)`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md)
  [`surv_optimize_allocation()`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md)
  : Optimize sequencing allocation across strata
- [`surv_plot_allocation()`](https://cuiweig.github.io/survinger/reference/surv_plot_allocation.md)
  : Plot allocation plan
- [`surv_plot_sequencing_rates()`](https://cuiweig.github.io/survinger/reference/surv_plot_sequencing_rates.md)
  : Plot sequencing rate inequality across strata
- [`surv_power_curve()`](https://cuiweig.github.io/survinger/reference/surv_power_curve.md)
  [`plot(`*`<surv_power_curve>`*`)`](https://cuiweig.github.io/survinger/reference/surv_power_curve.md)
  : Compute power curve for detection across prevalence range
- [`surv_prevalence_by()`](https://cuiweig.github.io/survinger/reference/surv_prevalence_by.md)
  : Estimate prevalence by subgroup
- [`surv_quality()`](https://cuiweig.github.io/survinger/reference/surv_quality.md)
  : Compute surveillance quality metrics
- [`surv_report()`](https://cuiweig.github.io/survinger/reference/surv_report.md)
  : Generate a comprehensive surveillance system report
- [`surv_reporting_probability()`](https://cuiweig.github.io/survinger/reference/surv_reporting_probability.md)
  : Compute cumulative reporting probability
- [`surv_required_sequences()`](https://cuiweig.github.io/survinger/reference/surv_required_sequences.md)
  : Required sequences for target detection probability
- [`surv_sensitivity()`](https://cuiweig.github.io/survinger/reference/surv_sensitivity.md)
  : Sensitivity analysis across methods
- [`surv_set_weights()`](https://cuiweig.github.io/survinger/reference/surv_set_weights.md)
  : Override design weights with custom values
- [`surv_simulate()`](https://cuiweig.github.io/survinger/reference/surv_simulate.md)
  : Simulate genomic surveillance data
- [`surv_table()`](https://cuiweig.github.io/survinger/reference/surv_table.md)
  : Format prevalence results for knitr tables
- [`surv_update_rates()`](https://cuiweig.github.io/survinger/reference/surv_update_rates.md)
  : Update sequencing rates in a surveillance design
- [`theme_survinger()`](https://cuiweig.github.io/survinger/reference/theme_survinger.md)
  : Publication-quality ggplot2 theme
- [`tidy(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  [`tidy(`*`<surv_nowcast>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  [`tidy(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  [`tidy(`*`<surv_allocation>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  [`tidy(`*`<surv_delay_fit>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  : Extract tidy estimates from survinger objects

## Other

- [`glance(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/glance.surv.md)
  [`glance(`*`<surv_delay_fit>`*`)`](https://cuiweig.github.io/survinger/reference/glance.surv.md)
  [`glance(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/glance.surv.md)
  : One-row summary of survinger model
- [`plot(`*`<surv_design>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_allocation>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_delay_fit>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_nowcast>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  [`plot(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/plot.surv.md)
  : Plot methods for survinger objects
- [`sarscov2_surveillance`](https://cuiweig.github.io/survinger/reference/sarscov2_surveillance.md)
  : Example SARS-CoV-2 genomic surveillance data
- [`print(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/surv_adjusted_prevalence.md)
  [`as.data.frame(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/surv_adjusted_prevalence.md)
  [`surv_adjusted_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_adjusted_prevalence.md)
  : Combined design-weighted and delay-adjusted prevalence
- [`surv_bind()`](https://cuiweig.github.io/survinger/reference/surv_bind.md)
  : Combine multiple prevalence estimates
- [`surv_compare_allocations()`](https://cuiweig.github.io/survinger/reference/surv_compare_allocations.md)
  : Compare multiple allocation strategies
- [`surv_compare_estimates()`](https://cuiweig.github.io/survinger/reference/surv_compare_estimates.md)
  : Compare weighted vs naive prevalence estimates
- [`print(`*`<surv_design>`*`)`](https://cuiweig.github.io/survinger/reference/surv_design.md)
  [`summary(`*`<surv_design>`*`)`](https://cuiweig.github.io/survinger/reference/surv_design.md)
  [`print(`*`<summary.surv_design>`*`)`](https://cuiweig.github.io/survinger/reference/surv_design.md)
  [`surv_design()`](https://cuiweig.github.io/survinger/reference/surv_design.md)
  : Create a genomic surveillance design object
- [`surv_design_effect()`](https://cuiweig.github.io/survinger/reference/surv_design_effect.md)
  : Compute design effect over time
- [`surv_detection_probability()`](https://cuiweig.github.io/survinger/reference/surv_detection_probability.md)
  : Variant detection probability under current design
- [`surv_estimate()`](https://cuiweig.github.io/survinger/reference/surv_estimate.md)
  : Pipe-friendly surveillance analysis
- [`print(`*`<surv_delay_fit>`*`)`](https://cuiweig.github.io/survinger/reference/surv_estimate_delay.md)
  [`surv_estimate_delay()`](https://cuiweig.github.io/survinger/reference/surv_estimate_delay.md)
  : Estimate reporting delay distribution
- [`surv_filter()`](https://cuiweig.github.io/survinger/reference/surv_filter.md)
  : Subset a surveillance design by filter criteria
- [`print(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md)
  [`as.data.frame(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md)
  [`surv_lineage_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_lineage_prevalence.md)
  : Estimate lineage prevalence with design weights
- [`surv_naive_prevalence()`](https://cuiweig.github.io/survinger/reference/surv_naive_prevalence.md)
  : Compute naive (unweighted) lineage prevalence
- [`print(`*`<surv_nowcast>`*`)`](https://cuiweig.github.io/survinger/reference/surv_nowcast_lineage.md)
  [`as.data.frame(`*`<surv_nowcast>`*`)`](https://cuiweig.github.io/survinger/reference/surv_nowcast_lineage.md)
  [`surv_nowcast_lineage()`](https://cuiweig.github.io/survinger/reference/surv_nowcast_lineage.md)
  : Nowcast lineage counts correcting for reporting delays
- [`print(`*`<surv_allocation>`*`)`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md)
  [`as.data.frame(`*`<surv_allocation>`*`)`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md)
  [`surv_optimize_allocation()`](https://cuiweig.github.io/survinger/reference/surv_optimize_allocation.md)
  : Optimize sequencing allocation across strata
- [`surv_plot_allocation()`](https://cuiweig.github.io/survinger/reference/surv_plot_allocation.md)
  : Plot allocation plan
- [`surv_plot_sequencing_rates()`](https://cuiweig.github.io/survinger/reference/surv_plot_sequencing_rates.md)
  : Plot sequencing rate inequality across strata
- [`surv_power_curve()`](https://cuiweig.github.io/survinger/reference/surv_power_curve.md)
  [`plot(`*`<surv_power_curve>`*`)`](https://cuiweig.github.io/survinger/reference/surv_power_curve.md)
  : Compute power curve for detection across prevalence range
- [`surv_prevalence_by()`](https://cuiweig.github.io/survinger/reference/surv_prevalence_by.md)
  : Estimate prevalence by subgroup
- [`surv_quality()`](https://cuiweig.github.io/survinger/reference/surv_quality.md)
  : Compute surveillance quality metrics
- [`surv_report()`](https://cuiweig.github.io/survinger/reference/surv_report.md)
  : Generate a comprehensive surveillance system report
- [`surv_reporting_probability()`](https://cuiweig.github.io/survinger/reference/surv_reporting_probability.md)
  : Compute cumulative reporting probability
- [`surv_required_sequences()`](https://cuiweig.github.io/survinger/reference/surv_required_sequences.md)
  : Required sequences for target detection probability
- [`surv_sensitivity()`](https://cuiweig.github.io/survinger/reference/surv_sensitivity.md)
  : Sensitivity analysis across methods
- [`surv_set_weights()`](https://cuiweig.github.io/survinger/reference/surv_set_weights.md)
  : Override design weights with custom values
- [`surv_simulate()`](https://cuiweig.github.io/survinger/reference/surv_simulate.md)
  : Simulate genomic surveillance data
- [`surv_table()`](https://cuiweig.github.io/survinger/reference/surv_table.md)
  : Format prevalence results for knitr tables
- [`surv_update_rates()`](https://cuiweig.github.io/survinger/reference/surv_update_rates.md)
  : Update sequencing rates in a surveillance design
- [`theme_survinger()`](https://cuiweig.github.io/survinger/reference/theme_survinger.md)
  : Publication-quality ggplot2 theme
- [`tidy(`*`<surv_prevalence>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  [`tidy(`*`<surv_nowcast>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  [`tidy(`*`<surv_adjusted>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  [`tidy(`*`<surv_allocation>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  [`tidy(`*`<surv_delay_fit>`*`)`](https://cuiweig.github.io/survinger/reference/tidy.surv.md)
  : Extract tidy estimates from survinger objects
