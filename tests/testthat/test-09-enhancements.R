test_that("tidy.surv_prevalence returns tibble", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 1)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  prev <- surv_lineage_prevalence(d, "BA.2.86")
  td <- tidy(prev)
  expect_s3_class(td, "tbl_df")
  expect_true("prevalence" %in% names(td))
})

test_that("glance.surv_prevalence returns one-row tibble", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 2)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  prev <- surv_lineage_prevalence(d, "BA.2.86")
  gl <- glance(prev)
  expect_equal(nrow(gl), 1)
  expect_true("mean_prevalence" %in% names(gl))
})

test_that("tidy.surv_allocation returns tibble", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 3)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  a <- surv_optimize_allocation(d, "min_mse", 200)
  td <- tidy(a)
  expect_s3_class(td, "tbl_df")
  expect_true("n_allocated" %in% names(td))
})

test_that("as.data.frame works for surv_prevalence", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 4)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  prev <- surv_lineage_prevalence(d, "BA.2.86")
  df <- as.data.frame(prev)
  expect_s3_class(df, "data.frame")
})

test_that("surv_filter subsets design correctly", {
  sim <- surv_simulate(n_regions = 5, n_weeks = 8, seed = 5)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  d_sub <- surv_filter(d, region %in% c("Region_A", "Region_B"))
  expect_s3_class(d_sub, "surv_design")
  expect_true(d_sub$n_obs < d$n_obs)
  expect_true(all(d_sub$data$region %in% c("Region_A", "Region_B")))
})

test_that("surv_report runs without error", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 6)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  result <- surv_report(d)
  expect_type(result, "list")
  expect_true("gini" %in% names(result))
  expect_true("detection_prob" %in% names(result))
})

test_that("lineage suggestion warning fires for misspelled lineage", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 7)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  expect_warning(
    surv_lineage_prevalence(d, "BA.2.8"),
    "not found"
  )
})

test_that("empty lineage filter returns empty nowcast", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 8)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  delay <- surv_estimate_delay(d)
  expect_warning(
    nc <- surv_nowcast_lineage(d, delay, "NONEXISTENT"),
    "No sequences"
  )
  expect_s3_class(nc, "surv_nowcast")
  expect_equal(nrow(nc$estimates), 0)
})

test_that("glance.surv_delay_fit works", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 9)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  fit <- surv_estimate_delay(d)
  gl <- glance(fit)
  expect_equal(nrow(gl), 1)
  expect_true("mean_delay" %in% names(gl))
})
