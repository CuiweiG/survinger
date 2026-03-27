test_that("surv_design creates valid object from data frame rates", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 4, seed = 1)
  design <- surv_design(
    data = sim$sequences, strata = ~ region,
    sequencing_rate = sim$population[c("region", "seq_rate")],
    population = sim$population
  )
  expect_s3_class(design, "surv_design")
  expect_equal(design$n_strata, 3)
  expect_true(design$n_obs > 0)
})

test_that("surv_design creates valid object from named vector rates", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 4, seed = 2)
  rates <- setNames(sim$population$seq_rate, sim$population$region)
  design <- surv_design(
    data = sim$sequences, strata = ~ region,
    sequencing_rate = rates, population = sim$population
  )
  expect_s3_class(design, "surv_design")
})

test_that("weights equal inverse of sequencing rate", {
  sim <- surv_simulate(n_regions = 4, n_weeks = 6, seed = 3)
  design <- surv_design(
    data = sim$sequences, strata = ~ region,
    sequencing_rate = sim$population[c("region", "seq_rate")],
    population = sim$population
  )
  expected <- 1 / design$strata_info$seq_rate
  expect_equal(design$weights$weight, expected, tolerance = 1e-10)
})

test_that("epiweek column is added", {
  sim <- surv_simulate(n_regions = 2, n_weeks = 3, seed = 4)
  design <- surv_design(
    data = sim$sequences, strata = ~ region,
    sequencing_rate = sim$population[c("region", "seq_rate")],
    population = sim$population
  )
  expect_true("epiweek" %in% names(design$data))
})

test_that("source_type is stored when provided", {
  sim <- surv_simulate(n_regions = 2, n_weeks = 3, seed = 5)
  design <- surv_design(
    data = sim$sequences, strata = ~ region,
    sequencing_rate = sim$population[c("region", "seq_rate")],
    population = sim$population, source_type = "source_type"
  )
  expect_equal(design$col_source_type, "source_type")
})

test_that("surv_update_rates recalculates weights", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 4, seed = 10)
  design <- surv_design(
    data = sim$sequences, strata = ~ region,
    sequencing_rate = sim$population[c("region", "seq_rate")],
    population = sim$population
  )
  old_wt <- design$weights$weight
  new_r <- sim$population[c("region", "seq_rate")]
  new_r$seq_rate <- new_r$seq_rate * 2
  updated <- surv_update_rates(design, new_r)
  expect_true(all(updated$weights$weight != old_wt))
})

test_that("surv_set_weights overrides", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 4, seed = 11)
  design <- surv_design(
    data = sim$sequences, strata = ~ region,
    sequencing_rate = sim$population[c("region", "seq_rate")],
    population = sim$population
  )
  updated <- surv_set_weights(design, rep(1.0, 3))
  expect_true(all(updated$weights$weight == 1.0))
})

test_that("surv_design rejects missing columns", {
  sim <- surv_simulate(n_regions = 2, n_weeks = 2, seed = 20)
  expect_error(
    surv_design(
      data = sim$sequences, strata = ~ nonexistent,
      sequencing_rate = sim$population[c("region", "seq_rate")],
      population = sim$population
    ), "not found"
  )
})

test_that("surv_design rejects non-formula strata", {
  sim <- surv_simulate(n_regions = 2, n_weeks = 2, seed = 21)
  expect_error(surv_design(
    data = sim$sequences, strata = "region",
    sequencing_rate = sim$population[c("region", "seq_rate")],
    population = sim$population
  ))
})

test_that("print and summary do not error", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 4, seed = 30)
  design <- surv_design(
    data = sim$sequences, strata = ~ region,
    sequencing_rate = sim$population[c("region", "seq_rate")],
    population = sim$population, source_type = "source_type"
  )
  expect_no_error(print(design))
  s <- summary(design)
  expect_s3_class(s, "summary.surv_design")
  expect_no_error(print(s))
})

test_that("surv_design with formula sequencing_rate works", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 4, seed = 40)
  pop <- sim$population
  design <- surv_design(
    data = sim$sequences, strata = ~ region,
    sequencing_rate = ~ n_sequenced / n_positive,
    population = pop
  )
  expect_s3_class(design, "surv_design")
  expect_true(all(design$strata_info$seq_rate > 0))
})
