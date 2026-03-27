# ============================================================
# Design-weighted lineage prevalence estimation
# ============================================================

#' Estimate lineage prevalence with design weights
#'
#' Estimates the prevalence of a specified pathogen lineage over time,
#' correcting for unequal sequencing rates across strata.
#'
#' @param design A `surv_design` object.
#' @param lineage Character. Target lineage name.
#' @param time Character. Time aggregation: `"epiweek"`, `"month"`,
#'   `"date"`, or a column name. Default `"epiweek"`.
#' @param method Character. Estimation method: `"hajek"` (default),
#'   `"horvitz_thompson"`, or `"poststratified"`.
#' @param conf_level Numeric. Confidence level. Default 0.95.
#' @param min_obs Integer. Minimum observations per time period.
#'   Default 5.
#'
#' @return A `surv_prevalence` object.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' prev <- surv_lineage_prevalence(d, "BA.2.86")
#' print(prev)
#'
#' @export
surv_lineage_prevalence <- function(design,
                                    lineage,
                                    time = "epiweek",
                                    method = c("hajek", "horvitz_thompson",
                                               "poststratified"),
                                    conf_level = 0.95,
                                    min_obs = 5L) {
  .assert_surv_design(design)
  checkmate::assert_string(lineage)
  method <- match.arg(method)
  checkmate::assert_number(conf_level, lower = 0.5, upper = 0.999)
  checkmate::assert_count(min_obs, positive = TRUE)

  dat <- .map_weights_to_obs(design)
  resolved <- .resolve_time_column(dat, time, design$col_date_collected)
  dat <- resolved$data
  time_col <- resolved$col

  dat$.is_target <- as.integer(dat[[design$col_lineage]] == lineage)
  z_alpha <- stats::qnorm(1 - (1 - conf_level) / 2)

  time_vals <- sort(unique(dat[[time_col]]))

  estimates <- purrr::map_dfr(time_vals, function(t) {
    d <- dat[dat[[time_col]] == t, , drop = FALSE]
    if (nrow(d) < min_obs) {
      return(tibble::tibble(
        time = t, lineage = lineage, n_obs = nrow(d),
        n_target = sum(d$.is_target), prevalence = NA_real_,
        se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        effective_n = NA_real_, flag = "insufficient_obs"
      ))
    }

    result <- switch(method,
      "horvitz_thompson" = .ht_estimator(d, design$strata_vars),
      "hajek"            = .hajek_estimator(d, design$strata_vars),
      "poststratified"   = .ps_estimator(d, design$strata_vars, design$strata_info)
    )

    tibble::tibble(
      time = t, lineage = lineage, n_obs = nrow(d),
      n_target = sum(d$.is_target),
      prevalence = result$estimate, se = result$se,
      ci_lower = max(0, result$estimate - z_alpha * result$se),
      ci_upper = min(1, result$estimate + z_alpha * result$se),
      effective_n = result$effective_n, flag = "ok"
    )
  })

  new_surv_prevalence(estimates, design, method, lineage, conf_level, time)
}


#' Compute naive (unweighted) lineage prevalence
#'
#' Simple proportion without design correction. Useful as baseline.
#'
#' @param design A `surv_design` object.
#' @param lineage Character. Target lineage.
#' @param time Character. Time aggregation. Default `"epiweek"`.
#' @param conf_level Numeric. Default 0.95.
#'
#' @return A `surv_prevalence` object with `method = "naive"`.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' naive <- surv_naive_prevalence(d, "BA.2.86")
#'
#' @export
surv_naive_prevalence <- function(design, lineage, time = "epiweek",
                                  conf_level = 0.95) {
  .assert_surv_design(design)
  checkmate::assert_string(lineage)

  dat <- design$data
  resolved <- .resolve_time_column(dat, time, design$col_date_collected)
  dat <- resolved$data
  time_col <- resolved$col
  lin_col <- design$col_lineage

  z_alpha <- stats::qnorm(1 - (1 - conf_level) / 2)
  time_vals <- sort(unique(dat[[time_col]]))

  estimates <- purrr::map_dfr(time_vals, function(t) {
    d <- dat[dat[[time_col]] == t, , drop = FALSE]
    n <- nrow(d)
    n_t <- sum(d[[lin_col]] == lineage)
    p <- n_t / max(n, 1L)
    se_val <- sqrt(p * (1 - p) / max(n, 1L))
    tibble::tibble(
      time = t, lineage = lineage, n_obs = n, n_target = n_t,
      prevalence = p, se = se_val,
      ci_lower = max(0, p - z_alpha * se_val),
      ci_upper = min(1, p + z_alpha * se_val),
      effective_n = as.double(n), flag = "ok"
    )
  })

  new_surv_prevalence(estimates, design, "naive", lineage, conf_level, time)
}


# ---- Internal estimators ----

#' @keywords internal
.ht_estimator <- function(d, strata_vars) {
  w <- d$weight
  y <- d$.is_target
  estimate <- sum(w * y) / sum(w)
  se <- .stratified_var_ht(d, strata_vars, y, w, estimate)
  eff_n <- .kish_eff_n(w)
  list(estimate = estimate, se = se, effective_n = eff_n)
}

#' @keywords internal
.hajek_estimator <- function(d, strata_vars) {
  w <- d$weight
  y <- d$.is_target
  denom <- sum(w)
  estimate <- sum(w * y) / denom
  residuals <- w * (y - estimate) / denom
  se <- .stratified_var_lin(d, strata_vars, residuals)
  eff_n <- .kish_eff_n(w)
  list(estimate = estimate, se = se, effective_n = eff_n)
}

#' @keywords internal
.ps_estimator <- function(d, strata_vars, strata_info) {
  stratum_prev <- d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(strata_vars))) |>
    dplyr::summarise(
      n_h = dplyr::n(), y_h = sum(.data$.is_target),
      p_h = .data$y_h / .data$n_h, .groups = "drop"
    )

  pop_col <- if ("pop_total" %in% names(strata_info)) "pop_total"
             else if ("n_positive" %in% names(strata_info)) "n_positive"
             else NULL

  if (!is.null(pop_col)) {
    stratum_prev <- dplyr::left_join(
      stratum_prev,
      strata_info[c(strata_vars, pop_col)],
      by = strata_vars
    )
    total_pop <- sum(strata_info[[pop_col]], na.rm = TRUE)
    stratum_prev$pi_h <- stratum_prev[[pop_col]] / total_pop
  } else {
    stratum_prev$pi_h <- 1 / nrow(stratum_prev)
  }

  estimate <- sum(stratum_prev$pi_h * stratum_prev$p_h)
  var_est <- sum(
    stratum_prev$pi_h^2 *
    stratum_prev$p_h * (1 - stratum_prev$p_h) /
    pmax(stratum_prev$n_h - 1, 1)
  )
  se <- sqrt(var_est)
  eff_n <- 1 / sum(stratum_prev$pi_h^2 / stratum_prev$n_h)
  list(estimate = estimate, se = se, effective_n = eff_n)
}


# ---- Variance helpers ----

#' @keywords internal
.stratified_var_ht <- function(d, strata_vars, y, w, estimate) {
  d$.y <- y
  d$.w <- w
  vc <- d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(strata_vars))) |>
    dplyr::summarise(
      n_h = dplyr::n(),
      f_h = 1 / mean(.data$.w),
      s2_h = if (dplyr::n() > 1L) stats::var(.data$.w * .data$.y) else 0,
      .groups = "drop"
    ) |>
    dplyr::mutate(v_h = .data$n_h * (1 - .data$f_h) / pmax(.data$n_h - 1, 1) * .data$s2_h)
  sqrt(sum(vc$v_h)) / sum(w)
}

#' @keywords internal
.stratified_var_lin <- function(d, strata_vars, residuals) {
  d$.resid <- residuals
  ve <- d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(strata_vars))) |>
    dplyr::summarise(
      n_h = dplyr::n(),
      v_h = .data$n_h / pmax(.data$n_h - 1, 1) *
            sum((.data$.resid - mean(.data$.resid))^2),
      .groups = "drop"
    )
  sqrt(sum(ve$v_h))
}
