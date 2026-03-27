# ============================================================
# Reporting delay estimation and nowcasting
# ============================================================

#' Estimate reporting delay distribution
#'
#' Fits a parametric distribution to the delay between sample collection
#' and sequence reporting, accounting for right-truncation.
#'
#' @param design A `surv_design` object with both collection and report dates.
#' @param distribution Character: `"negbin"` (default), `"poisson"`,
#'   `"lognormal"`, or `"nonparametric"`.
#' @param strata Optional one-sided formula for delay stratification.
#' @param max_delay Integer. Maximum plausible delay in days. Default 60.
#' @param ref_date Date. Reference date for right-truncation. Default
#'   is max report date in data.
#'
#' @return A `surv_delay_fit` object.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 12, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' fit <- surv_estimate_delay(d)
#' print(fit)
#'
#' @export
surv_estimate_delay <- function(design,
                                distribution = c("negbin", "poisson",
                                                 "lognormal", "nonparametric"),
                                strata = NULL,
                                max_delay = 60L,
                                ref_date = NULL) {

  .assert_surv_design(design)
  distribution <- match.arg(distribution)
  checkmate::assert_count(max_delay, positive = TRUE)

  if (is.null(design$col_date_reported)) {
    cli::cli_abort("Design has no reporting date column. Cannot estimate delays.")
  }

  dat <- design$data
  d_col <- design$col_date_collected
  r_col <- design$col_date_reported

  dat$.delay <- as.integer(dat[[r_col]] - dat[[d_col]])
  dat <- dat[dat$.delay >= 0 & dat$.delay <= max_delay, , drop = FALSE]

  if (nrow(dat) == 0L) {
    cli::cli_abort("No valid delay observations after filtering.")
  }

  if (is.null(ref_date)) ref_date <- max(dat[[r_col]], na.rm = TRUE)
  ref_date <- as.Date(ref_date)
  dat$.trunc_time <- as.integer(ref_date - dat[[d_col]])

  if (!is.null(strata)) {
    strata_vars_delay <- all.vars(strata)
    .assert_columns_exist(dat, strata_vars_delay)
  } else {
    strata_vars_delay <- NULL
  }

  if (is.null(strata_vars_delay)) {
    params <- .fit_delay_single(dat$.delay, dat$.trunc_time,
                                distribution, max_delay)
    parameters <- tibble::tibble(stratum = "all", distribution = distribution,
                                 !!!params)
  } else {
    groups <- split(dat, dat[strata_vars_delay])
    parameters <- purrr::map_dfr(names(groups), function(g) {
      d_sub <- groups[[g]]
      p <- .fit_delay_single(d_sub$.delay, d_sub$.trunc_time,
                             distribution, max_delay)
      tibble::tibble(stratum = g, distribution = distribution, !!!p)
    })
  }

  new_surv_delay_fit(
    distribution = distribution,
    parameters = parameters,
    strata = strata_vars_delay,
    data_summary = list(n = nrow(dat), mean_delay = mean(dat$.delay),
                        median_delay = stats::median(dat$.delay)),
    diagnostics = list(ref_date = ref_date, max_delay = max_delay,
                       prop_truncated = mean(dat$.delay >= dat$.trunc_time))
  )
}


#' Compute cumulative reporting probability
#'
#' @param delay_fit A `surv_delay_fit` object.
#' @param delta Integer vector. Days since collection.
#' @param stratum Character or `NULL`. Default `NULL`.
#'
#' @return Numeric vector of probabilities.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 12, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' fit <- surv_estimate_delay(d)
#' surv_reporting_probability(fit, delta = c(7, 14, 21))
#'
#' @export
surv_reporting_probability <- function(delay_fit, delta, stratum = NULL) {
  .assert_surv_class(delay_fit, "surv_delay_fit")
  checkmate::assert_integerish(delta, lower = 0)

  params <- .get_delay_params(delay_fit, stratum)
  dist <- delay_fit$distribution

  if (dist == "nonparametric") {
    cdf <- params$cdf
    idx <- pmin(as.integer(delta) + 1L, length(cdf))
    return(cdf[idx])
  }

  switch(dist,
    "negbin"    = stats::pnbinom(delta, mu = params$mu, size = params$size),
    "poisson"   = stats::ppois(delta, lambda = params$lambda),
    "lognormal" = stats::plnorm(delta + 1, meanlog = params$meanlog,
                                sdlog = params$sdlog)
  )
}


#' Nowcast lineage counts correcting for reporting delays
#'
#' @param design A `surv_design` object.
#' @param delay_fit A `surv_delay_fit` object.
#' @param lineage Character or `NULL`. Target lineage.
#' @param time Character. Default `"epiweek"`.
#' @param horizon Integer. Recent periods to nowcast. Default 4.
#' @param ref_date Date or `NULL`. Default max report date.
#' @param method Character: `"direct"` (default) or `"em"`.
#'
#' @return A `surv_nowcast` object.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 12, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' fit <- surv_estimate_delay(d)
#' nc <- surv_nowcast_lineage(d, fit, "BA.2.86")
#' print(nc)
#'
#' @export
surv_nowcast_lineage <- function(design, delay_fit, lineage = NULL,
                                 time = "epiweek", horizon = 4L,
                                 ref_date = NULL,
                                 method = c("direct", "em")) {
  .assert_surv_design(design)
  .assert_surv_class(delay_fit, "surv_delay_fit")
  method <- match.arg(method)
  checkmate::assert_count(horizon, positive = TRUE)

  dat <- design$data
  d_col <- design$col_date_collected
  r_col <- design$col_date_reported
  lin_col <- design$col_lineage

  if (is.null(ref_date)) {
    ref_date <- if (!is.null(r_col)) max(dat[[r_col]], na.rm = TRUE) else Sys.Date()
  }
  ref_date <- as.Date(ref_date)

  resolved <- .resolve_time_column(dat, time, d_col)
  dat <- resolved$data
  time_col <- resolved$col

  if (!is.null(lineage)) dat <- dat[dat[[lin_col]] == lineage, , drop = FALSE]

  time_summary <- dat |>
    dplyr::mutate(.time = .data[[time_col]]) |>
    dplyr::group_by(.data$.time) |>
    dplyr::summarise(
      n_observed = dplyr::n(),
      median_collection = stats::median(.data[[d_col]]),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$.time) |>
    dplyr::mutate(days_ago = as.integer(ref_date - .data$median_collection))

  time_summary$report_prob <- surv_reporting_probability(
    delay_fit, pmax(time_summary$days_ago, 0L)
  )
  time_summary$report_prob <- pmax(time_summary$report_prob, 0.01)

  all_times <- sort(unique(time_summary$.time))
  nowcast_times <- utils::tail(all_times, horizon)

  if (method == "direct") {
    estimates <- time_summary |>
      dplyr::mutate(
        n_estimated = .data$n_observed / .data$report_prob,
        se = sqrt(.data$n_estimated * (1 - .data$report_prob) /
                    .data$report_prob),
        ci_lower = pmax(0, .data$n_estimated - stats::qnorm(0.975) * .data$se),
        ci_upper = .data$n_estimated + stats::qnorm(0.975) * .data$se,
        is_nowcast = .data$.time %in% nowcast_times
      ) |>
      dplyr::rename(time = ".time")
  } else {
    estimates <- .nowcast_em(time_summary, nowcast_times)
  }

  new_surv_nowcast(estimates, delay_fit, horizon, method, lineage)
}


# ---- Internal delay fitting ----

#' @keywords internal
.fit_delay_single <- function(delays, trunc_times, distribution, max_delay) {
  if (distribution == "nonparametric") {
    return(.fit_delay_nonparametric(delays, trunc_times, max_delay))
  }

  negloglik <- function(par) {
    if (distribution == "negbin") {
      mu <- exp(par[1]); size <- exp(par[2])
      log_f <- stats::dnbinom(delays, mu = mu, size = size, log = TRUE)
      log_F <- stats::pnbinom(trunc_times, mu = mu, size = size, log.p = TRUE)
    } else if (distribution == "poisson") {
      lambda <- exp(par[1])
      log_f <- stats::dpois(delays, lambda = lambda, log = TRUE)
      log_F <- stats::ppois(trunc_times, lambda = lambda, log.p = TRUE)
    } else {
      meanlog <- par[1]; sdlog <- exp(par[2])
      log_f <- log(pmax(stats::plnorm(delays + 1, meanlog, sdlog) -
                          stats::plnorm(delays, meanlog, sdlog), 1e-300))
      log_F <- stats::plnorm(trunc_times + 1, meanlog, sdlog, log.p = TRUE)
    }
    -sum(log_f - log_F)
  }

  init <- switch(distribution,
    "negbin"    = c(log(mean(delays) + 0.1), log(1)),
    "poisson"   = log(mean(delays) + 0.1),
    "lognormal" = c(log(mean(delays) + 1), log(0.5))
  )

  fit <- stats::optim(init, negloglik, method = "Nelder-Mead",
                      control = list(maxit = 5000))

  switch(distribution,
    "negbin"    = list(mu = exp(fit$par[1]), size = exp(fit$par[2]),
                       converged = fit$convergence == 0),
    "poisson"   = list(lambda = exp(fit$par[1]),
                       converged = fit$convergence == 0),
    "lognormal" = list(meanlog = fit$par[1], sdlog = exp(fit$par[2]),
                       converged = fit$convergence == 0)
  )
}

#' @keywords internal
.fit_delay_nonparametric <- function(delays, trunc_times, max_delay) {
  pmf <- rep(1 / (max_delay + 1), max_delay + 1)
  for (iter in seq_len(50)) {
    cdf_vals <- cumsum(pmf)
    inc_prob <- cdf_vals[pmin(trunc_times + 1L, max_delay + 1L)]
    inc_prob <- pmax(inc_prob, 1e-10)
    weights <- 1 / inc_prob
    new_pmf <- rep(0, max_delay + 1)
    for (i in seq_along(delays)) {
      idx <- delays[i] + 1L
      if (idx <= max_delay + 1L) new_pmf[idx] <- new_pmf[idx] + weights[i]
    }
    new_pmf <- new_pmf / sum(new_pmf)
    if (max(abs(new_pmf - pmf)) < 1e-8) break
    pmf <- new_pmf
  }
  list(pmf = pmf, cdf = cumsum(pmf), converged = TRUE)
}

#' @keywords internal
.get_delay_params <- function(delay_fit, stratum = NULL) {
  params <- delay_fit$parameters
  if (is.null(stratum)) {
    if ("all" %in% params$stratum) {
      return(as.list(params[params$stratum == "all", ][1, ]))
    }
    return(as.list(params[1, ]))
  }
  row <- params[params$stratum == stratum, ]
  if (nrow(row) == 0L) cli::cli_abort("Stratum not found in delay fit.")
  as.list(row[1, ])
}

#' @keywords internal
.nowcast_em <- function(time_summary, nowcast_times,
                        max_iter = 200, tol = 1e-6) {
  n_obs <- time_summary$n_observed
  rp <- time_summary$report_prob
  n_est <- n_obs / rp

  for (iter in seq_len(max_iter)) {
    n_new <- n_obs + n_est * (1 - rp)
    if (max(abs(n_new - n_est) / pmax(n_est, 1)) < tol) break
    n_est <- n_new
  }

  se <- sqrt(n_est * (1 - rp) / rp^2)

  time_summary |>
    dplyr::mutate(
      n_estimated = n_est, se = se,
      ci_lower = pmax(0, n_est - stats::qnorm(0.975) * se),
      ci_upper = n_est + stats::qnorm(0.975) * se,
      is_nowcast = .data$.time %in% nowcast_times
    ) |>
    dplyr::rename(time = ".time")
}
