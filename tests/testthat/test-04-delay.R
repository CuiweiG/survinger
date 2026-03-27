test_that("delay estimation recovers known negbin parameters", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 20,
                       delay_params = list(mu = 10, size = 3), seed = 99)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  fit <- surv_estimate_delay(d, distribution = "negbin")
  expect_s3_class(fit, "surv_delay_fit")
  expect_equal(fit$parameters$mu[1], 10, tolerance = 3)
  expect_equal(fit$parameters$size[1], 3, tolerance = 2)
})

test_that("reporting probability is monotonically increasing", {
  sim <- surv_simulate(n_regions = 2, n_weeks = 10, seed = 50)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  fit <- surv_estimate_delay(d)
  probs <- surv_reporting_probability(fit, delta = 0:60)
  expect_true(all(diff(probs) >= -1e-10))
})

test_that("reporting probability approaches 1", {
  sim <- surv_simulate(n_regions = 2, n_weeks = 10, seed = 51)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  fit <- surv_estimate_delay(d)
  expect_true(surv_reporting_probability(fit, delta = 120L) > 0.99)
})

test_that("nowcast inflates recent counts", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 12, seed = 60)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  fit <- surv_estimate_delay(d)
  nc <- surv_nowcast_lineage(d, fit, "BA.2.86")
  expect_s3_class(nc, "surv_nowcast")
  recent <- nc$estimates[nc$estimates$is_nowcast, ]
  if (nrow(recent) > 0) {
    expect_true(all(recent$n_estimated >= recent$n_observed - 0.1))
  }
})

test_that("direct and em give similar results", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 12, seed = 61)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  fit <- surv_estimate_delay(d)
  nc_d <- surv_nowcast_lineage(d, fit, "BA.2.86", method = "direct")
  nc_e <- surv_nowcast_lineage(d, fit, "BA.2.86", method = "em")
  diff_val <- abs(nc_d$estimates$n_estimated - nc_e$estimates$n_estimated)
  expect_true(all(diff_val / pmax(nc_d$estimates$n_estimated, 1) < 0.5))
})

test_that("print.surv_delay_fit works", {
  sim <- surv_simulate(n_regions = 2, n_weeks = 8, seed = 62)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  fit <- surv_estimate_delay(d)
  expect_no_error(print(fit))
})
