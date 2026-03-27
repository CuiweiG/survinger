test_that("equal allocation divides evenly", {
  sim <- surv_simulate(n_regions = 4, n_weeks = 4, seed = 1)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  n <- .equal_allocation(4, 100)
  expect_equal(sum(n), 100)
  expect_equal(length(n), 4)
})

test_that("total allocation equals capacity for all objectives", {
  sim <- surv_simulate(n_regions = 4, n_weeks = 8, seed = 2)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  for (obj in c("min_mse", "max_detection", "min_imbalance")) {
    a <- surv_optimize_allocation(d, obj, total_capacity = 200,
                                  target_prevalence = 0.05)
    expect_equal(sum(a$allocation$n_allocated), 200)
  }
})

test_that("min_per_stratum is respected", {
  sim <- surv_simulate(n_regions = 5, n_weeks = 8, seed = 3)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  a <- surv_optimize_allocation(d, "min_mse", total_capacity = 100,
                                min_per_stratum = 10)
  expect_true(all(a$allocation$n_allocated >= 10))
})

test_that("allocation is integer-valued", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 6, seed = 4)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  a <- surv_optimize_allocation(d, "min_mse", total_capacity = 500)
  expect_true(all(a$allocation$n_allocated == as.integer(a$allocation$n_allocated)))
})

test_that("compare_allocations returns correct structure", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 6, seed = 5)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  comp <- surv_compare_allocations(d, total_capacity = 200)
  expect_s3_class(comp, "tbl_df")
  expect_true("strategy" %in% names(comp))
  expect_true(nrow(comp) == 5)
})

test_that("print.surv_allocation works", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 4, seed = 6)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  a <- surv_optimize_allocation(d, "min_mse", total_capacity = 100)
  expect_no_error(print(a))
})
