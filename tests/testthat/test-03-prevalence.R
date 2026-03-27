test_that("naive prevalence matches simple proportion", {
  sim <- surv_simulate(n_regions = 2, n_weeks = 4, seed = 10)
  design <- surv_design(sim$sequences, ~ region,
                        sim$population[c("region", "seq_rate")], sim$population)
  naive <- surv_naive_prevalence(design, "BA.2.86", time = "epiweek")
  first_week <- naive$estimates$time[1]
  d <- sim$sequences[sim$sequences$epiweek == first_week, ]
  expected <- sum(d$lineage == "BA.2.86") / nrow(d)
  expect_equal(naive$estimates$prevalence[1], expected, tolerance = 1e-10)
})

test_that("weighted differs from naive under unequal rates", {
  sim <- surv_simulate(n_regions = 5, n_weeks = 12,
                       sequencing_rates = c(0.5, 0.05, 0.02, 0.01, 0.3),
                       seed = 42)
  design <- surv_design(sim$sequences, ~ region,
                        sim$population[c("region", "seq_rate")], sim$population)
  weighted <- surv_lineage_prevalence(design, "BA.2.86", method = "hajek")
  naive <- surv_naive_prevalence(design, "BA.2.86")
  diffs <- abs(weighted$estimates$prevalence - naive$estimates$prevalence)
  expect_true(any(diffs > 0.005, na.rm = TRUE))
})

test_that("all three methods return valid objects", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 20)
  design <- surv_design(sim$sequences, ~ region,
                        sim$population[c("region", "seq_rate")], sim$population)
  for (m in c("hajek", "horvitz_thompson", "poststratified")) {
    prev <- surv_lineage_prevalence(design, "BA.2.86", method = m)
    expect_s3_class(prev, "surv_prevalence")
    valid <- !is.na(prev$estimates$prevalence)
    expect_true(all(prev$estimates$prevalence[valid] >= 0))
    expect_true(all(prev$estimates$prevalence[valid] <= 1))
  }
})

test_that("effective_n < n_obs when weights vary", {
  sim <- surv_simulate(n_regions = 4, n_weeks = 8,
                       sequencing_rates = c(0.5, 0.05, 0.01, 0.005), seed = 30)
  design <- surv_design(sim$sequences, ~ region,
                        sim$population[c("region", "seq_rate")], sim$population)
  prev <- surv_lineage_prevalence(design, "BA.2.86")
  valid <- !is.na(prev$estimates$effective_n)
  expect_true(all(prev$estimates$effective_n[valid] <= prev$estimates$n_obs[valid]))
})

test_that("min_obs filtering works", {
  sim <- surv_simulate(n_regions = 2, n_weeks = 2, seed = 40)
  design <- surv_design(sim$sequences, ~ region,
                        sim$population[c("region", "seq_rate")], sim$population)
  prev <- surv_lineage_prevalence(design, "BA.2.86", min_obs = 99999L)
  expect_true(all(is.na(prev$estimates$prevalence)))
  expect_true(all(prev$estimates$flag == "insufficient_obs"))
})

test_that("CI contains point estimate", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 50)
  design <- surv_design(sim$sequences, ~ region,
                        sim$population[c("region", "seq_rate")], sim$population)
  prev <- surv_lineage_prevalence(design, "BA.2.86")
  valid <- prev$estimates$flag == "ok"
  est <- prev$estimates[valid, ]
  expect_true(all(est$prevalence >= est$ci_lower))
  expect_true(all(est$prevalence <= est$ci_upper))
})

test_that("hajek matches survey::svymean for pooled estimate", {
  skip_if_not_installed("survey")
  sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 100)
  design <- surv_design(sim$sequences, ~ region,
                        sim$population[c("region", "seq_rate")], sim$population)
  dat <- design$data
  dat$.target <- as.integer(dat$lineage == "BA.2.86")
  wt_info <- design$strata_info[c("region", "seq_rate")]
  dat <- dplyr::left_join(dat, wt_info, by = "region", relationship = "many-to-one")
  dat$wt <- 1 / dat$seq_rate
  survinger_est <- sum(dat$wt * dat$.target) / sum(dat$wt)
  svy <- survey::svydesign(ids = ~1, strata = ~region, weights = ~wt, data = dat)
  survey_est <- as.numeric(coef(survey::svymean(~.target, svy)))
  expect_equal(survinger_est, survey_est, tolerance = 0.001)
})

test_that("poststratified estimator works", {
  sim <- surv_simulate(n_regions = 3, n_weeks = 8, seed = 90)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")], sim$population)
  res <- surv_lineage_prevalence(d, "BA.5", method = "poststratified")
  expect_s3_class(res, "surv_prevalence")
})
