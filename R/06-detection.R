# ============================================================
# Detection probability
# ============================================================

#' Variant detection probability under current design
#'
#' @param design A `surv_design` object.
#' @param true_prevalence Numeric in (0,1).
#' @param delay_fit Optional `surv_delay_fit`.
#' @param n_periods Integer. Accumulation periods. Default 1.
#' @param detection_threshold Integer. Min detections. Default 1.
#'
#' @return A list with `overall`, `cumulative`, `by_stratum`, `parameters`.
#'
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' surv_detection_probability(d, 0.01)
#'
#' @export
surv_detection_probability <- function(design, true_prevalence,
                                       delay_fit = NULL,
                                       n_periods = 1L,
                                       detection_threshold = 1L) {
  .assert_surv_design(design)
  checkmate::assert_number(true_prevalence, lower = 0, upper = 1)
  checkmate::assert_count(n_periods, positive = TRUE)
  checkmate::assert_count(detection_threshold, positive = TRUE)

  info <- design$strata_info
  n_weeks <- length(unique(design$data[["epiweek"]]))
  n_seq <- info$n_sequenced_obs / max(1L, n_weeks)

  p_detect <- 1 - (1 - true_prevalence)^n_seq

  if (!is.null(delay_fit)) {
    avg_rp <- surv_reporting_probability(delay_fit, delta = 7L)
    p_detect <- p_detect * avg_rp
  }

  lambda <- -sum(log(1 - pmin(p_detect, 0.999)))
  p_per <- if (detection_threshold == 1L) {
    1 - exp(-lambda)
  } else {
    stats::ppois(detection_threshold - 1L, lambda, lower.tail = FALSE)
  }

  p_cum <- numeric(n_periods)
  p_no <- 1
  for (k in seq_len(n_periods)) {
    p_no <- p_no * (1 - p_per)
    p_cum[k] <- 1 - p_no
  }

  list(
    overall = p_per,
    cumulative = tibble::tibble(period = seq_len(n_periods), p_detect = p_cum),
    by_stratum = tibble::tibble(
      stratum_id = info$stratum_id, n_seq_per_period = n_seq,
      p_detect = p_detect
    ),
    parameters = list(true_prevalence = true_prevalence,
                      n_periods = n_periods,
                      detection_threshold = detection_threshold,
                      delay_adjusted = !is.null(delay_fit))
  )
}


#' Required sequences for target detection probability
#'
#' @param true_prevalence Numeric.
#' @param target_probability Numeric. Default 0.95.
#' @param n_periods Integer. Default 1.
#' @param detection_threshold Integer. Default 1.
#'
#' @return Integer.
#'
#' @examples
#' surv_required_sequences(0.01)
#' surv_required_sequences(0.05, target_probability = 0.99)
#'
#' @export
surv_required_sequences <- function(true_prevalence,
                                    target_probability = 0.95,
                                    n_periods = 1L,
                                    detection_threshold = 1L) {
  checkmate::assert_number(true_prevalence, lower = 1e-6, upper = 1)
  checkmate::assert_number(target_probability, lower = 0, upper = 1)

  if (detection_threshold == 1L && n_periods == 1L) {
    return(as.integer(ceiling(
      log(1 - target_probability) / log(1 - true_prevalence)
    )))
  }

  for (n in seq_len(100000L)) {
    lambda <- -n * log(1 - true_prevalence)
    p_per <- stats::ppois(detection_threshold - 1L, lambda, lower.tail = FALSE)
    p_cum <- 1 - (1 - p_per)^n_periods
    if (p_cum >= target_probability) return(as.integer(n))
  }
  NA_integer_
}
