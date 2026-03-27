test_that("detection probability increases with prevalence", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 5)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  p1 <- surv_detection_probability(d, 0.01)
  p2 <- surv_detection_probability(d, 0.05)
  p3 <- surv_detection_probability(d, 0.10)
  expect_true(p1$overall < p2$overall)
  expect_true(p2$overall < p3$overall)
})

test_that("required sequences decreases with prevalence", {
  n1 <- surv_required_sequences(0.01)
  n2 <- surv_required_sequences(0.05)
  n3 <- surv_required_sequences(0.10)
  expect_true(n1 > n2)
  expect_true(n2 > n3)
})

test_that("required_sequences gives correct known result", {
  n <- surv_required_sequences(0.01, target_probability = 0.95)
  expect_equal(n, 299L)
})

test_that("detection probability with delay adjustment", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 6)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  delay <- surv_estimate_delay(d)
  det <- surv_detection_probability(d, 0.05, delay_fit = delay)
  expect_true(det$overall > 0 && det$overall < 1)
})

test_that("multi-period detection accumulates", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 7)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  n1 <- surv_detection_probability(d, 0.01, n_periods = 1)
  n4 <- surv_detection_probability(d, 0.01, n_periods = 4)
  expect_true(n4$cumulative$p_detect[4] > n1$overall)
})
