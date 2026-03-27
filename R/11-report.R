# ============================================================
# Comprehensive surveillance report
# ============================================================

#' Generate a comprehensive surveillance system report
#'
#' Produces a summary of the current surveillance design's
#' strengths, weaknesses, and recommendations.
#'
#' @param design A `surv_design` object.
#' @param target_lineage Character or `NULL`. Lineage to focus on.
#'   If `NULL`, uses the most common non-"Other" lineage.
#' @param target_prevalence Numeric. Assumed prevalence for detection
#'   calculations. Default 0.01.
#'
#' @return Invisibly returns a named list of computed metrics including
#'   n_obs, n_strata, rate_range, gini, effective_n, detection_prob,
#'   and mean_bias.
#'
#' @seealso [surv_design()], [surv_lineage_prevalence()],
#'   [surv_detection_probability()], [surv_optimize_allocation()]
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' surv_report(d)
#'
#' @export
surv_report <- function(design,
                        target_lineage = NULL,
                        target_prevalence = 0.01) {
  .assert_surv_design(design)

  info <- design$strata_info
  dat <- design$data

  # Auto-select target lineage
  if (is.null(target_lineage)) {
    lin_tab <- sort(table(dat[[design$col_lineage]]), decreasing = TRUE)
    non_other <- names(lin_tab)[!grepl("^other$", names(lin_tab),
                                       ignore.case = TRUE)]
    target_lineage <- if (length(non_other) > 0) non_other[1] else names(lin_tab)[1]
  }

  n_obs <- design$n_obs
  n_strata <- design$n_strata
  rate_range <- range(info$seq_rate)
  rate_ratio <- rate_range[2] / rate_range[1]
  gini <- .gini_coefficient(info$seq_rate)

  weighted_dat <- .map_weights_to_obs(design)
  eff_n <- .kish_eff_n(weighted_dat$weight)
  eff_ratio <- eff_n / n_obs

  det <- surv_detection_probability(design, target_prevalence)
  req <- surv_required_sequences(target_prevalence)

  weighted <- surv_lineage_prevalence(design, target_lineage, method = "hajek")
  naive <- surv_naive_prevalence(design, target_lineage)
  bias <- mean(abs(weighted$estimates$prevalence - naive$estimates$prevalence),
               na.rm = TRUE)

  cli::cli_h1("Surveillance System Report")

  cli::cli_h2("Design Structure")
  n_fmt <- formatC(n_obs, big.mark = ",")
  cli::cli_text("Total sequences: {.strong {n_fmt}}")
  cli::cli_text("Strata: {.strong {n_strata}}")
  dr <- range(dat[[design$col_date_collected]], na.rm = TRUE)
  cli::cli_text("Period: {dr[1]} to {dr[2]}")

  cli::cli_h2("Sequencing Inequality")
  lo <- round(rate_range[1] * 100, 2)
  hi <- round(rate_range[2] * 100, 2)
  rr <- round(rate_ratio, 1)
  gi <- round(gini, 3)
  cli::cli_text("Rate range: {lo}% to {hi}%")
  cli::cli_text("Rate ratio: {.strong {rr}}x")
  cli::cli_text("Gini coefficient: {.strong {gi}}")
  if (gini > 0.4) {
    cli::cli_alert_warning("High inequality. Design weighting is essential.")
  } else if (gini > 0.2) {
    cli::cli_alert_info("Moderate inequality. Weighting recommended.")
  } else {
    cli::cli_alert_success("Low inequality. Naive estimates may be acceptable.")
  }

  cli::cli_h2("Estimation Impact")
  en <- round(eff_n)
  er <- round(eff_ratio * 100, 1)
  bp <- round(bias * 100, 2)
  cli::cli_text("Effective sample size: {.strong {en}} ({er}% of total)")
  cli::cli_text("Mean |weighted - naive| bias: {.strong {bp}} percentage points")

  cli::cli_h2("Detection Power")
  dp <- round(det$overall * 100, 1)
  tp <- target_prevalence * 100
  cli::cli_text("Target: detect {target_lineage} at {tp}% prevalence")
  cli::cli_text("Weekly detection probability: {.strong {dp}}%")
  cli::cli_text("Sequences needed for 95% detection: {.strong {req}}")

  invisible(list(
    n_obs = n_obs, n_strata = n_strata,
    rate_range = rate_range, rate_ratio = rate_ratio, gini = gini,
    effective_n = eff_n, eff_ratio = eff_ratio,
    detection_prob = det$overall, required_sequences = req,
    mean_bias = bias
  ))
}

#' @keywords internal
.gini_coefficient <- function(x) {
  x <- sort(x)
  n <- length(x)
  if (n <= 1L || sum(x) == 0) return(0)
  i <- seq_len(n)
  2 * sum(i * x) / (n * sum(x)) - (n + 1) / n
}
