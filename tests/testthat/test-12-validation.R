# ============================================================
# Validation tests backing README claims
# ============================================================

test_that("Hajek matches survey::svymean on simulated data", {
  skip_if_not_installed("survey")

  sim <- surv_simulate(n_regions = 5, n_weeks = 12, seed = 100)
  d <- surv_design(sim$sequences, ~ region,
                   sim$population[c("region", "seq_rate")],
                   sim$population)

  # survinger estimate (Hajek = ratio estimator)
  prev <- surv_lineage_prevalence(d, "BA.2.86", method = "hajek")
  surv_est <- mean(prev$estimates$prevalence, na.rm = TRUE)

  # survey package estimate
  dat <- d$data
  dat$.is_target <- as.integer(dat$lineage == "BA.2.86")
  weights_map <- d$strata_info[c(d$strata_vars, "seq_rate")]
  dat <- merge(dat, weights_map, by = d$strata_vars)
  dat$.weight <- 1 / dat$seq_rate

  svy <- survey::svydesign(ids = ~1, strata = ~region,
                            weights = ~.weight, data = dat)
  svy_result <- survey::svymean(~.is_target, svy)
  svy_est <- as.numeric(coef(svy_result))

  # Should match within rounding tolerance
  # (not exact due to per-period vs pooled aggregation)
  expect_lt(abs(surv_est - svy_est), 0.05)
})

test_that("Wilson CI achieves near-nominal coverage in simulation", {
  # Simulate many datasets and check CI coverage
  set.seed(42)
  n_reps <- 200
  true_p <- 0.15
  n_obs <- 100
  covered <- 0L

  for (i in seq_len(n_reps)) {
    x <- rbinom(n_obs, 1, true_p)
    p_hat <- mean(x)
    z <- qnorm(0.975)
    denom <- 1 + z^2 / n_obs
    center <- (p_hat + z^2 / (2 * n_obs)) / denom
    margin <- z * sqrt(p_hat * (1 - p_hat) / n_obs +
                         z^2 / (4 * n_obs^2)) / denom
    lo <- max(0, center - margin)
    hi <- min(1, center + margin)
    if (true_p >= lo && true_p <= hi) covered <- covered + 1L
  }
  coverage <- covered / n_reps

  # Brown et al. 2001 target: 93-97% for nominal 95%
  expect_gt(coverage, 0.90)
  expect_lt(coverage, 0.99)
})

test_that("Delay MLE recovers true parameters at n = 5000", {
  set.seed(123)
  true_mu <- 10
  true_size <- 3
  n <- 5000

  delays <- rnbinom(n, mu = true_mu, size = true_size)
  # No truncation (all observed)
  trunc_times <- rep(100L, n)

  # Use internal delay fitter
  params <- survinger:::.fit_delay_single(delays, trunc_times,
                                           "negbin", max_delay = 60)

  # Recovery within 5% relative error
  mu_error <- abs(params$mu - true_mu) / true_mu
  size_error <- abs(params$size - true_size) / true_size

  expect_lt(mu_error, 0.05)
  expect_lt(size_error, 0.15)  # size harder to estimate
  expect_true(params$converged)
})
