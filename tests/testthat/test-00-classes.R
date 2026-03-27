test_that("new_surv_prevalence creates correct class", {
  obj <- new_surv_prevalence(
    estimates = tibble::tibble(time = "2024-W01", prevalence = 0.3),
    design = NULL, method = "hajek", lineage = "BA.2.86",
    conf_level = 0.95, time_unit = "epiweek"
  )
  expect_s3_class(obj, "surv_prevalence")
  expect_equal(obj$method, "hajek")
  expect_equal(obj$lineage, "BA.2.86")
})

test_that("new_surv_allocation creates correct class", {
  obj <- new_surv_allocation(
    allocation = tibble::tibble(stratum = "A", n_allocated = 100),
    objective = "min_mse", total_capacity = 100,
    constraints = list(budget = NULL), diagnostics = list(converged = TRUE)
  )
  expect_s3_class(obj, "surv_allocation")
  expect_equal(obj$objective, "min_mse")
})

test_that("new_surv_delay_fit creates correct class", {
  obj <- new_surv_delay_fit(
    distribution = "negbin",
    parameters = tibble::tibble(stratum = "all", mu = 10, size = 3),
    strata = NULL,
    data_summary = list(n = 100, mean_delay = 10, median_delay = 8),
    diagnostics = list()
  )
  expect_s3_class(obj, "surv_delay_fit")
})

test_that("new_surv_nowcast creates correct class", {
  obj <- new_surv_nowcast(
    estimates = tibble::tibble(time = "2024-W01", n_observed = 50),
    delay_fit = NULL, truncation_window = 4,
    method = "direct", lineage = "BA.2.86"
  )
  expect_s3_class(obj, "surv_nowcast")
})

test_that("new_surv_adjusted creates correct class", {
  obj <- new_surv_adjusted(
    estimates = tibble::tibble(time = "2024-W01", prevalence = 0.3),
    design_component = NULL, delay_component = NULL,
    method = "hajek+direct"
  )
  expect_s3_class(obj, "surv_adjusted")
})

test_that("print methods do not error", {
  prev <- new_surv_prevalence(
    estimates = tibble::tibble(time = "2024-W01", prevalence = 0.3,
                               se = 0.05, ci_lower = 0.2, ci_upper = 0.4),
    design = NULL, method = "hajek", lineage = "BA.2.86",
    conf_level = 0.95, time_unit = "epiweek"
  )
  expect_no_error(print(prev))

  alloc <- new_surv_allocation(
    allocation = tibble::tibble(stratum = "A", n_allocated = 100),
    objective = "min_mse", total_capacity = 100,
    constraints = list(), diagnostics = list()
  )
  expect_no_error(print(alloc))

  delay <- new_surv_delay_fit(
    distribution = "negbin",
    parameters = tibble::tibble(stratum = "all", mu = 10, size = 3),
    strata = NULL,
    data_summary = list(n = 100, mean_delay = 10, median_delay = 8),
    diagnostics = list()
  )
  expect_no_error(print(delay))
})
