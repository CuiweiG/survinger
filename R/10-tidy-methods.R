# ============================================================
# Tidy methods (broom-style) for surv_* objects
# ============================================================

#' Extract tidy estimates from survinger objects
#'
#' Converts survinger result objects into tidy tibbles suitable for
#' further analysis with dplyr, ggplot2, or other tidyverse tools.
#'
#' @param x A survinger result object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A tibble.
#'
#' @name tidy.surv
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' prev <- surv_lineage_prevalence(d, "BA.2.86")
#' tidy(prev)
NULL

#' @rdname tidy.surv
#' @export
tidy.surv_prevalence <- function(x, ...) {
  x$estimates
}

#' @rdname tidy.surv
#' @export
tidy.surv_nowcast <- function(x, ...) {
  x$estimates
}

#' @rdname tidy.surv
#' @export
tidy.surv_adjusted <- function(x, ...) {
  x$estimates
}

#' @rdname tidy.surv
#' @export
tidy.surv_allocation <- function(x, ...) {
  x$allocation
}

#' @rdname tidy.surv
#' @export
tidy.surv_delay_fit <- function(x, ...) {
  x$parameters
}


#' One-row summary of survinger model
#'
#' @param x A survinger result object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A single-row tibble with key summary statistics.
#'
#' @name glance.surv
#' @examples
#' sim <- surv_simulate(n_regions = 3, n_weeks = 10, seed = 1)
#' d <- surv_design(sim$sequences, ~ region,
#'                  sim$population[c("region", "seq_rate")], sim$population)
#' prev <- surv_lineage_prevalence(d, "BA.2.86")
#' glance(prev)
NULL

#' @rdname glance.surv
#' @export
glance.surv_prevalence <- function(x, ...) {
  est <- x$estimates[!is.na(x$estimates$prevalence), , drop = FALSE]
  tibble::tibble(
    lineage = x$lineage,
    method = x$method,
    n_periods = nrow(est),
    mean_prevalence = mean(est$prevalence),
    median_prevalence = stats::median(est$prevalence),
    mean_se = mean(est$se, na.rm = TRUE),
    mean_effective_n = mean(est$effective_n, na.rm = TRUE),
    conf_level = x$conf_level
  )
}

#' @rdname glance.surv
#' @export
glance.surv_delay_fit <- function(x, ...) {
  tibble::tibble(
    distribution = x$distribution,
    n_observations = x$data_summary$n,
    mean_delay = x$data_summary$mean_delay,
    median_delay = x$data_summary$median_delay
  )
}

#' @rdname glance.surv
#' @export
glance.surv_adjusted <- function(x, ...) {
  est <- x$estimates[!is.na(x$estimates$prevalence), , drop = FALSE]
  tibble::tibble(
    method = x$method,
    n_periods = nrow(est),
    mean_prevalence = mean(est$prevalence),
    mean_se = mean(est$se, na.rm = TRUE),
    mean_report_prob = mean(est$mean_report_prob, na.rm = TRUE)
  )
}
