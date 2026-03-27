test_that("adjusted prevalence runs end-to-end", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 12, seed = 77)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  delay <- surv_estimate_delay(d)
  adj <- surv_adjusted_prevalence(d, delay, "BA.2.86")
  expect_s3_class(adj, "surv_adjusted")
  expect_true(nrow(adj$estimates) > 0)
  vals <- adj$estimates$prevalence[!is.na(adj$estimates$prevalence)]
  expect_true(all(vals >= 0 & vals <= 1))
})

test_that("adjusted has both components", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 78)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  delay <- surv_estimate_delay(d)
  adj <- surv_adjusted_prevalence(d, delay, "BA.2.86")
  expect_s3_class(adj$design_component, "surv_prevalence")
  expect_s3_class(adj$delay_component, "surv_nowcast")
})

test_that("print.surv_adjusted works", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 79)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  delay <- surv_estimate_delay(d)
  adj <- surv_adjusted_prevalence(d, delay, "BA.2.86")
  expect_no_error(print(adj))
})
