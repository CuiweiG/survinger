# ============================================================
# S3 plot methods for all surv_* classes
# ============================================================

#' @export
plot.surv_design <- function(x, ...) {
  surv_plot_sequencing_rates(x)
}

#' @export
plot.surv_allocation <- function(x, ...) {
  surv_plot_allocation(x)
}

#' @export
plot.surv_prevalence <- function(x, ...) {
  est <- x$estimates[!is.na(x$estimates$prevalence), , drop = FALSE]
  ggplot2::ggplot(est, ggplot2::aes(x = .data$time, y = .data$prevalence, group = 1)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      fill = "#2E75B6", alpha = 0.2
    ) +
    ggplot2::geom_line(color = "#2E75B6", linewidth = 0.8) +
    ggplot2::geom_point(color = "#2E75B6", size = 1.5) +
    ggplot2::scale_y_continuous(
      labels = function(v) paste0(round(v * 100, 1), "%"),
      limits = c(0, NA)
    ) +
    ggplot2::labs(title = paste0(x$lineage, " prevalence (", x$method, ")"),
                  x = "Time", y = "Prevalence") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' @export
plot.surv_delay_fit <- function(x, ...) {
  params <- x$parameters[1, ]
  d_seq <- 0:60
  pmf_vals <- switch(x$distribution,
    "negbin"    = stats::dnbinom(d_seq, mu = params$mu, size = params$size),
    "poisson"   = stats::dpois(d_seq, lambda = params$lambda),
    "lognormal" = stats::plnorm(d_seq + 1, params$meanlog, params$sdlog) -
                  stats::plnorm(d_seq, params$meanlog, params$sdlog),
    "nonparametric" = {
      pmf <- if (is.list(params$pmf)) params$pmf[[1]] else params$pmf
      pmf_out <- rep(0, 61)
      pmf_out[seq_len(min(length(pmf), 61))] <- pmf[seq_len(min(length(pmf), 61))]
      pmf_out
    }
  )
  if (length(pmf_vals) < 61) {
    pmf_vals <- c(pmf_vals, rep(0, 61 - length(pmf_vals)))
  }

  df <- tibble::tibble(delay = d_seq, probability = pmf_vals[seq_len(61)])
  ggplot2::ggplot(df, ggplot2::aes(x = .data$delay, y = .data$probability)) +
    ggplot2::geom_col(fill = "#D4652F", alpha = 0.7) +
    ggplot2::labs(title = paste0("Reporting Delay (", x$distribution, ")"),
                  x = "Delay (days)", y = "Probability") +
    ggplot2::theme_minimal(base_size = 12)
}

#' @export
plot.surv_nowcast <- function(x, ...) {
  est <- x$estimates
  ggplot2::ggplot(est, ggplot2::aes(x = .data$time, group = 1)) +
    ggplot2::geom_col(ggplot2::aes(y = .data$n_observed),
                      fill = "grey70", alpha = 0.5) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      fill = "#D4652F", alpha = 0.2
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$n_estimated),
                       color = "#D4652F", linewidth = 0.8) +
    ggplot2::labs(
      title = paste0("Nowcast",
                     if (!is.null(x$lineage)) paste0(": ", x$lineage) else ""),
      x = "Time", y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' @export
plot.surv_adjusted <- function(x, ...) {
  est <- x$estimates[!is.na(x$estimates$prevalence), , drop = FALSE]
  lin_label <- if (nrow(est) > 0 && "lineage" %in% names(est)) est$lineage[1] else ""
  ggplot2::ggplot(est, ggplot2::aes(x = .data$time, group = 1)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      fill = "#5B8C5A", alpha = 0.2
    ) +
    ggplot2::geom_line(ggplot2::aes(y = .data$prevalence),
                       color = "#5B8C5A", linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = .data$prevalence),
                        color = "#5B8C5A", size = 1.5) +
    ggplot2::scale_y_continuous(
      labels = function(v) paste0(round(v * 100, 1), "%"),
      limits = c(0, NA)
    ) +
    ggplot2::labs(title = paste0(lin_label, " \u2014 Design + Delay Adjusted"),
                  x = "Time", y = "Prevalence") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
