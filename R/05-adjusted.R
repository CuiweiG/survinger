# ============================================================
# Combined design-weighted + delay-adjusted inference
# ============================================================

#' Combined design-weighted and delay-adjusted prevalence
#'
#' Simultaneously corrects for unequal sequencing rates and
#' right-truncation from reporting delays.
#'
#' @param design A `surv_design` object.
#' @param delay_fit A `surv_delay_fit` object.
#' @param lineage Character. Target lineage.
#' @param time Character. Default `"epiweek"`.
#' @param prevalence_method Character. Default `"hajek"`.
#' @param nowcast_method Character. Default `"direct"`.
#' @param conf_level Numeric. Default 0.95.
#' @param bootstrap_n Integer. 0 for delta method, >0 for bootstrap.
#'   Default 0.
#'
#' @return A `surv_adjusted` object.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 12, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' delay <- surv_estimate_delay(d)
#' adj <- surv_adjusted_prevalence(d, delay, "BA.2.86")
#' print(adj)
#'
#' @export
surv_adjusted_prevalence <- function(design, delay_fit, lineage,
                                     time = "epiweek",
                                     prevalence_method = "hajek",
                                     nowcast_method = "direct",
                                     conf_level = 0.95,
                                     bootstrap_n = 0L) {
  .assert_surv_design(design)
  .assert_surv_class(delay_fit, "surv_delay_fit")
  checkmate::assert_string(lineage)
  checkmate::assert_count(bootstrap_n)

  dat <- .map_weights_to_obs(design)
  d_col <- design$col_date_collected
  lin_col <- design$col_lineage
  r_col <- design$col_date_reported
  sv <- design$strata_vars

  ref_date <- if (!is.null(r_col)) max(dat[[r_col]], na.rm = TRUE) else Sys.Date()
  ref_date <- as.Date(ref_date)

  resolved <- .resolve_time_column(dat, time, d_col)
  dat <- resolved$data
  time_col <- resolved$col

  z_alpha <- stats::qnorm(1 - (1 - conf_level) / 2)
  dat$.is_target <- as.integer(dat[[lin_col]] == lineage)
  dat$.time <- dat[[time_col]]

  # Per time-stratum cell: compute report probability
  cell <- dat |>
    dplyr::group_by(.data$.time, !!!rlang::syms(sv)) |>
    dplyr::summarise(
      n_obs = dplyr::n(), n_target = sum(.data$.is_target),
      mean_weight = mean(.data$weight),
      med_coll = stats::median(.data[[d_col]]),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      days_ago = as.integer(ref_date - .data$med_coll),
      report_prob = surv_reporting_probability(delay_fit, pmax(.data$days_ago, 0L)),
      report_prob = pmax(.data$report_prob, 0.01),
      n_obs_adj = .data$n_obs / .data$report_prob,
      n_target_adj = .data$n_target / .data$report_prob
    )

  # Aggregate by time with design weights
  estimates <- cell |>
    dplyr::group_by(.data$.time) |>
    dplyr::summarise(
      n_obs_raw = sum(.data$n_obs),
      n_obs_adjusted = sum(.data$n_obs_adj),
      prevalence = sum(.data$mean_weight * .data$n_target_adj) /
        sum(.data$mean_weight * .data$n_obs_adj),
      var_design = {
        p <- sum(.data$mean_weight * .data$n_target) /
          sum(.data$mean_weight * .data$n_obs)
        sum(.data$mean_weight^2 * .data$n_obs * p * (1 - p)) /
          sum(.data$mean_weight * .data$n_obs)^2
      },
      var_delay = sum(.data$mean_weight^2 * .data$n_target *
                        (1 - .data$report_prob) /
                        pmax(.data$report_prob^2, 1e-4)) /
        sum(.data$mean_weight * .data$n_obs_adj)^2,
      mean_report_prob = stats::weighted.mean(.data$report_prob, .data$n_obs),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      se = sqrt(.data$var_design + .data$var_delay),
      ci_lower = pmax(0, .data$prevalence - z_alpha * .data$se),
      ci_upper = pmin(1, .data$prevalence + z_alpha * .data$se),
      lineage = lineage
    ) |>
    dplyr::rename(time = ".time") |>
    dplyr::select("time", "lineage", "n_obs_raw", "n_obs_adjusted",
                  "prevalence", "se", "ci_lower", "ci_upper",
                  "mean_report_prob")

  design_comp <- surv_lineage_prevalence(design, lineage, time,
                                         prevalence_method, conf_level)
  delay_comp <- surv_nowcast_lineage(design, delay_fit, lineage, time,
                                     method = nowcast_method)

  new_surv_adjusted(estimates, design_comp, delay_comp,
                    paste0("design:", prevalence_method, "+delay:", nowcast_method))
}
