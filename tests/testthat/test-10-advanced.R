test_that("surv_prevalence_by produces per-group estimates", {
  sim <- surv_simulate(n_regions = 4, n_weeks = 8, seed = 1)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  by_region <- surv_prevalence_by(d, "BA.2.86", by = "region")
  expect_s3_class(by_region, "tbl_df")
  expect_true("group" %in% names(by_region))
  expect_true(length(unique(by_region$group)) >= 3)
})

test_that("surv_sensitivity returns all methods", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 2)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  sens <- surv_sensitivity(d, "BA.2.86")
  expect_s3_class(sens, "tbl_df")
  expect_true(all(c("naive", "hajek", "poststratified") %in% sens$method))
})

test_that("surv_sensitivity includes adjusted when delay_fit provided", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 3)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  delay <- surv_estimate_delay(d)
  sens <- surv_sensitivity(d, "BA.2.86", delay_fit = delay)
  expect_true("adjusted" %in% sens$method)
})

test_that("surv_power_curve returns correct structure", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 4)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  pc <- surv_power_curve(d)
  expect_s3_class(pc, "surv_power_curve")
  expect_true(nrow(pc$curve) > 0)
  expect_true(nrow(pc$thresholds) == 3)
  expect_true(all(pc$curve$detection >= 0 & pc$curve$detection <= 1))
})

test_that("surv_power_curve plot works", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 5)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  pc <- surv_power_curve(d)
  p <- plot(pc)
  expect_s3_class(p, "ggplot")
})

test_that("surv_power_curve print works", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 6)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  pc <- surv_power_curve(d)
  expect_no_error(print(pc))
})
