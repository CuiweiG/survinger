test_that("surv_compare_estimates returns ggplot", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 1)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  w <- surv_lineage_prevalence(d, "BA.2.86")
  n <- surv_naive_prevalence(d, "BA.2.86")
  p <- surv_compare_estimates(w, n)
  expect_s3_class(p, "ggplot")
})

test_that("surv_plot_sequencing_rates returns ggplot", {
  sim <- surv_simulate(n_regions = 4, n_weeks = 6, seed = 2)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  p <- surv_plot_sequencing_rates(d)
  expect_s3_class(p, "ggplot")
})

test_that("surv_plot_allocation returns ggplot", {
  sim <- surv_simulate(n_regions = 4, n_weeks = 6, seed = 3)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  a <- surv_optimize_allocation(d, "min_mse", total_capacity = 200)
  p <- surv_plot_allocation(a)
  expect_s3_class(p, "ggplot")
})

test_that("surv_design_effect returns tibble", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 4)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  w <- surv_lineage_prevalence(d, "BA.2.86")
  n <- surv_naive_prevalence(d, "BA.2.86")
  deff <- surv_design_effect(w, n)
  expect_s3_class(deff, "tbl_df")
  expect_true("deff" %in% names(deff))
  expect_true("bias_correction" %in% names(deff))
})

test_that("all plot.surv_* methods run without error", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 5)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)

  expect_s3_class(plot(d), "ggplot")

  a <- surv_optimize_allocation(d, "min_mse", 200)
  expect_s3_class(plot(a), "ggplot")

  w <- surv_lineage_prevalence(d, "BA.2.86")
  expect_s3_class(plot(w), "ggplot")

  delay <- surv_estimate_delay(d)
  expect_s3_class(plot(delay), "ggplot")

  nc <- surv_nowcast_lineage(d, delay, "BA.2.86")
  expect_s3_class(plot(nc), "ggplot")

  adj <- surv_adjusted_prevalence(d, delay, "BA.2.86")
  expect_s3_class(plot(adj), "ggplot")
})
