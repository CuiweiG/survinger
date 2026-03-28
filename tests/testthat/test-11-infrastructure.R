test_that("surv_estimate one-liner works", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 1)
  result <- surv_estimate(
    data = sim$sequences, strata = ~ region,
    sequencing_rate = sim$population[c("region", "seq_rate")],
    population = sim$population, lineage = "BA.2.86"
  )
  expect_s3_class(result, "surv_prevalence")
})

test_that("surv_estimate with delay correction", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 2)
  result <- surv_estimate(
    data = sim$sequences, strata = ~ region,
    sequencing_rate = sim$population[c("region", "seq_rate")],
    population = sim$population, lineage = "BA.2.86",
    correct_delay = TRUE
  )
  expect_s3_class(result, "surv_adjusted")
})

test_that("surv_bind combines multiple results", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 3)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  p1 <- surv_lineage_prevalence(d, "BA.5")
  p2 <- surv_lineage_prevalence(d, "XBB.1.5")
  combined <- surv_bind(p1, p2)
  expect_s3_class(combined, "tbl_df")
  expect_true("source" %in% names(combined))
  expect_true(length(unique(combined$source)) == 2)
})

test_that("surv_table formats prevalence", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 4)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  prev <- surv_lineage_prevalence(d, "BA.2.86")
  tbl <- surv_table(prev)
  expect_s3_class(tbl, "tbl_df")
  expect_true("ci" %in% names(tbl))
})

test_that("surv_table formats allocation", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 5)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  alloc <- surv_optimize_allocation(d, "min_mse", 200)
  tbl <- surv_table(alloc)
  expect_s3_class(tbl, "tbl_df")
})

test_that("surv_quality returns one-row tibble", {
  sim <- surv_simulate(n_regions = 4, n_weeks = 10, seed = 6)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  q <- surv_quality(d)
  expect_s3_class(q, "tbl_df")
  expect_equal(nrow(q), 1)
  expect_true(all(c("gini", "deff", "rate_ratio", "detection_at_1pct") %in% names(q)))
  expect_true(q$deff >= 1)
})
